package main

import (
	"crypto/tls"
	"crypto/x509"
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/Shopify/sarama"
	"github.com/krallistic/kazoo-go"
	"github.com/pkg/errors"
	"github.com/prometheus/client_golang/prometheus"
	"github.com/prometheus/client_golang/prometheus/promhttp"
	plog "github.com/prometheus/common/promlog"
	plogflag "github.com/prometheus/common/promlog/flag"
	"github.com/prometheus/common/version"
	"github.com/rcrowley/go-metrics"
	"gopkg.in/alecthomas/kingpin.v2"
	"k8s.io/klog/v2"
)

const (
	namespace = "fast5"
	clientID  = "fast5_exporter"
)

const (
	INFO  = 0
	DEBUG = 1
	TRACE = 2
)

var (
	totalSizeMetric           *prometheus.Desc
	totalAverageQualityMetric *prometheus.Desc
	numberOfReadsMetric       *prometheus.Desc
	nMetric                   *prometheus.Desc
	aMetric                   *prometheus.Desc
	cMetric                   *prometheus.Desc
	gMetric                   *prometheus.Desc
	tMetric                   *prometheus.Desc
	maxReadMetric             *prometheus.Desc
	averageQualityMetric      *prometheus.Desc
	rawDataLengthMetric       *prometheus.Desc
)

var filenameSizeMap = make(map[string]int64)

var numberOfReadsMap = make(map[string]int)
var numberOfNucleotidesMap = make(map[string]int)
var numberOfNucleotidesTotal = 0
var averageQualityTotal = float64(0.0)
var numberOfNMap = make(map[string]int)
var numberOfAMap = make(map[string]int)
var numberOfCMap = make(map[string]int)
var numberOfGMap = make(map[string]int)
var numberOfTMap = make(map[string]int)
var maxReadMap = make(map[string]int)
var averageQualityMap = make(map[string]float64)
var rawDataLengthMap = make(map[string]int)

type Exporter struct {
	client                  sarama.Client
	topicFilter             *regexp.Regexp
	groupFilter             *regexp.Regexp
	mu                      sync.Mutex
	useZooKeeperLag         bool
	zookeeperClient         *kazoo.Kazoo
	nextMetadataRefresh     time.Time
	metadataRefreshInterval time.Duration
	offsetShowAll           bool
	topicWorkers            int
	allowConcurrent         bool
	sgMutex                 sync.Mutex
	sgWaitCh                chan struct{}
	sgChans                 []chan<- prometheus.Metric
	consumerGroupFetchAll   bool
}

type exporterOpts struct {
	uri                      []string
	useSASL                  bool
	useSASLHandshake         bool
	saslUsername             string
	saslPassword             string
	saslMechanism            string
	saslDisablePAFXFast      bool
	useTLS                   bool
	tlsServerName            string
	tlsCAFile                string
	tlsCertFile              string
	tlsKeyFile               string
	serverUseTLS             bool
	serverMutualAuthEnabled  bool
	serverTlsCAFile          string
	serverTlsCertFile        string
	serverTlsKeyFile         string
	tlsInsecureSkipTLSVerify bool
	fast5Version             string
	useZooKeeperLag          bool
	uriZookeeper             []string
	labels                   string
	metadataRefreshInterval  string
	serviceName              string
	kerberosConfigPath       string
	realm                    string
	keyTabPath               string
	kerberosAuthType         string
	offsetShowAll            bool
	topicWorkers             int
	allowConcurrent          bool
	allowAutoTopicCreation   bool
	verbosityLogLevel        int
}

// CanReadCertAndKey returns true if the certificate and key files already exists,
// otherwise returns false. If lost one of cert and key, returns error.
func CanReadCertAndKey(certPath, keyPath string) (bool, error) {
	certReadable := canReadFile(certPath)
	keyReadable := canReadFile(keyPath)

	if certReadable == false && keyReadable == false {
		return false, nil
	}

	if certReadable == false {
		return false, fmt.Errorf("error reading %s, certificate and key must be supplied as a pair", certPath)
	}

	if keyReadable == false {
		return false, fmt.Errorf("error reading %s, certificate and key must be supplied as a pair", keyPath)
	}

	return true, nil
}

// If the file represented by path exists and
// readable, returns true otherwise returns false.
func canReadFile(path string) bool {
	f, err := os.Open(path)
	if err != nil {
		return false
	}

	defer f.Close()

	return true
}

// NewExporter returns an initialized Exporter.
func NewExporter(opts exporterOpts, topicFilter string, groupFilter string) (*Exporter, error) {
	var zookeeperClient *kazoo.Kazoo
	config := sarama.NewConfig()
	config.ClientID = clientID
	fast5Version, err := sarama.ParseKafkaVersion(opts.fast5Version)
	if err != nil {
		return nil, err
	}
	config.Version = fast5Version

	if opts.useSASL {
		// Convert to lowercase so that SHA512 and SHA256 is still valid
		opts.saslMechanism = strings.ToLower(opts.saslMechanism)
		switch opts.saslMechanism {
		case "scram-sha512":
			config.Net.SASL.SCRAMClientGeneratorFunc = func() sarama.SCRAMClient { return &XDGSCRAMClient{HashGeneratorFcn: SHA512} }
			config.Net.SASL.Mechanism = sarama.SASLMechanism(sarama.SASLTypeSCRAMSHA512)
		case "scram-sha256":
			config.Net.SASL.SCRAMClientGeneratorFunc = func() sarama.SCRAMClient { return &XDGSCRAMClient{HashGeneratorFcn: SHA256} }
			config.Net.SASL.Mechanism = sarama.SASLMechanism(sarama.SASLTypeSCRAMSHA256)
		case "gssapi":
			config.Net.SASL.Mechanism = sarama.SASLMechanism(sarama.SASLTypeGSSAPI)
			config.Net.SASL.GSSAPI.ServiceName = opts.serviceName
			config.Net.SASL.GSSAPI.KerberosConfigPath = opts.kerberosConfigPath
			config.Net.SASL.GSSAPI.Realm = opts.realm
			config.Net.SASL.GSSAPI.Username = opts.saslUsername
			if opts.kerberosAuthType == "keytabAuth" {
				config.Net.SASL.GSSAPI.AuthType = sarama.KRB5_KEYTAB_AUTH
				config.Net.SASL.GSSAPI.KeyTabPath = opts.keyTabPath
			} else {
				config.Net.SASL.GSSAPI.AuthType = sarama.KRB5_USER_AUTH
				config.Net.SASL.GSSAPI.Password = opts.saslPassword
			}
			if opts.saslDisablePAFXFast {
				config.Net.SASL.GSSAPI.DisablePAFXFAST = true
			}
		case "plain":
		default:
			return nil, fmt.Errorf(
				`invalid sasl mechanism "%s": can only be "scram-sha256", "scram-sha512", "gssapi" or "plain"`,
				opts.saslMechanism,
			)
		}

		config.Net.SASL.Enable = true
		config.Net.SASL.Handshake = opts.useSASLHandshake

		if opts.saslUsername != "" {
			config.Net.SASL.User = opts.saslUsername
		}

		if opts.saslPassword != "" {
			config.Net.SASL.Password = opts.saslPassword
		}
	}

	if opts.useTLS {
		config.Net.TLS.Enable = true

		config.Net.TLS.Config = &tls.Config{
			ServerName:         opts.tlsServerName,
			RootCAs:            x509.NewCertPool(),
			InsecureSkipVerify: opts.tlsInsecureSkipTLSVerify,
		}

		if opts.tlsCAFile != "" {
			if ca, err := ioutil.ReadFile(opts.tlsCAFile); err == nil {
				config.Net.TLS.Config.RootCAs.AppendCertsFromPEM(ca)
			} else {
				return nil, err
			}
		}

		canReadCertAndKey, err := CanReadCertAndKey(opts.tlsCertFile, opts.tlsKeyFile)
		if err != nil {
			return nil, errors.Wrap(err, "error reading cert and key")
		}
		if canReadCertAndKey {
			cert, err := tls.LoadX509KeyPair(opts.tlsCertFile, opts.tlsKeyFile)
			if err == nil {
				config.Net.TLS.Config.Certificates = []tls.Certificate{cert}
			} else {
				return nil, err
			}
		}
	}

	if opts.useZooKeeperLag {
		klog.V(DEBUG).Infoln("Using zookeeper lag, so connecting to zookeeper")
		zookeeperClient, err = kazoo.NewKazoo(opts.uriZookeeper, nil)
		if err != nil {
			return nil, errors.Wrap(err, "error connecting to zookeeper")
		}
	}

	interval, err := time.ParseDuration(opts.metadataRefreshInterval)
	if err != nil {
		return nil, errors.Wrap(err, "Cannot parse metadata refresh interval")
	}

	config.Metadata.RefreshFrequency = interval

	config.Metadata.AllowAutoTopicCreation = opts.allowAutoTopicCreation

	client, err := sarama.NewClient(opts.uri, config)

	if err != nil {
		return nil, errors.Wrap(err, "Error Init Fast5 Client")
	}

	klog.V(TRACE).Infoln("Done Init Clients")
	// Init our exporter.
	return &Exporter{
		client:                  client,
		topicFilter:             regexp.MustCompile(topicFilter),
		groupFilter:             regexp.MustCompile(groupFilter),
		useZooKeeperLag:         opts.useZooKeeperLag,
		zookeeperClient:         zookeeperClient,
		nextMetadataRefresh:     time.Now(),
		metadataRefreshInterval: interval,
		offsetShowAll:           opts.offsetShowAll,
		topicWorkers:            opts.topicWorkers,
		allowConcurrent:         opts.allowConcurrent,
		sgMutex:                 sync.Mutex{},
		sgWaitCh:                nil,
		sgChans:                 []chan<- prometheus.Metric{},
		consumerGroupFetchAll:   config.Version.IsAtLeast(sarama.V2_0_0_0),
	}, nil
}

//func (e *Exporter) fetchOffsetVersion() int16 {
//	version := e.client.Config().Version
//	if e.client.Config().Version.IsAtLeast(sarama.V2_0_0_0) {
//		return 4
//	} else if version.IsAtLeast(sarama.V0_10_2_0) {
//		return 2
//	} else if version.IsAtLeast(sarama.V0_8_2_2) {
//		return 1
//	}
//	return 0
//}

// Describe describes all the metrics ever exported by the Fast5 exporter. It
// implements prometheus.Collector.
func (e *Exporter) Describe(ch chan<- *prometheus.Desc) {
	ch <- totalSizeMetric
	ch <- numberOfReadsMetric
	ch <- nMetric
	ch <- aMetric
	ch <- cMetric
	ch <- gMetric
	ch <- tMetric
	ch <- maxReadMetric
	ch <- averageQualityMetric
	ch <- rawDataLengthMetric
}

func (e *Exporter) Collect(ch chan<- prometheus.Metric) {
	//	   	if e.allowConcurrent {
	e.collect(ch)
	return
	//	   	}
	// Locking to avoid race add
	/*e.sgMutex.Lock()
	  e.sgChans = append(e.sgChans, ch)
	  // Safe to compare length since we own the Lock

	  	if len(e.sgChans) == 1 {
	  		e.sgWaitCh = make(chan struct{})
	  		go e.collectChans(e.sgWaitCh)
	  	} else {

	  		klog.V(TRACE).Info("concurrent calls detected, waiting for first to finish")
	  	}

	  // Put in another variable to ensure not overwriting it in another Collect once we wait
	  waiter := e.sgWaitCh
	  e.sgMutex.Unlock()
	  // Released lock, we have insurance that our chan will be part of the collectChan slice
	  <-waiter*/
	// collectChan finished

}

func (e *Exporter) collectChans(quit chan struct{}) {
	original := make(chan prometheus.Metric)
	container := make([]prometheus.Metric, 0, 100)
	go func() {
		for metric := range original {
			container = append(container, metric)
		}
	}()
	e.collect(original)
	close(original)
	e.sgMutex.Lock()
	for _, ch := range e.sgChans {
		for _, metric := range container {
			ch <- metric
		}
	}
	// Reset the slice
	e.sgChans = e.sgChans[:0]
	close(quit)
	e.sgMutex.Unlock()
}

type statisticsData struct {
	FileSize      string `json:'fileSize'`
	RawDataLength string `json:'rawDataLength'`
	ChannelNumber string `json:'channelNumber'`
}

func (e *Exporter) collect(ch chan<- prometheus.Metric) {
	filepath.Walk("/tmp/fast5", func(path string, info os.FileInfo, err error) error {
		if filepath.Ext(path) == ".fast5" {
			if _, ok := filenameSizeMap[path]; ok {
				if filenameSizeMap[path] == info.Size() {
					//continue
				} else {
					filenameSizeMap[path] = info.Size()
					e.run(path)
				}

			} else {
				filenameSizeMap[path] = info.Size()
				e.run(path)
			}
		}

		if err != nil {
			fmt.Println("ERROR:", err)
		}
		return nil
	})

	var totalSize = 0
	for _, element := range filenameSizeMap {
		totalSize += int(element)
	}
	ch <- prometheus.MustNewConstMetric(
		totalSizeMetric, prometheus.GaugeValue, float64(totalSize), "address", "name",
	)
	/*ch <- prometheus.MustNewConstMetric(
		totalAverageQualityMetric, prometheus.GaugeValue, float64(averageQualityTotal), "address", "name",
	)

	for key, element := range numberOfReadsMap {
		ch <- prometheus.MustNewConstMetric(
			numberOfReadsMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range numberOfNMap {
		ch <- prometheus.MustNewConstMetric(
			nMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range numberOfAMap {
		ch <- prometheus.MustNewConstMetric(
			aMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range numberOfCMap {
		ch <- prometheus.MustNewConstMetric(
			cMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range numberOfGMap {
		ch <- prometheus.MustNewConstMetric(
			gMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range numberOfTMap {
		ch <- prometheus.MustNewConstMetric(
			tMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range maxReadMap {
		ch <- prometheus.MustNewConstMetric(
			maxReadMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
	for key, element := range averageQualityMap {
		ch <- prometheus.MustNewConstMetric(
			averageQualityMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}*/
	for key, element := range rawDataLengthMap {
		ch <- prometheus.MustNewConstMetric(
			rawDataLengthMetric, prometheus.GaugeValue, float64(element), "address", "name", key,
		)
	}
}

func (e *Exporter) run(path string) {
	var statistics []*statisticsData
	result, err := exec.Command("python3", "python/parse_fast5_file.py", "--path", path).Output()
	if err != nil {
		log.Fatal(err)
	}
	err = json.Unmarshal([]byte(result), &statistics)
	if err != nil {
		log.Fatal(err)
	}
	for _, statistic := range statistics {
		RawDataLength, _ := strconv.Atoi(statistic.RawDataLength)
		/*NumberOfReads, _ := strconv.Atoi(statistic.NumberOfReads)
		NumberOfN, _ := strconv.Atoi(statistic.NumberOfN)
		NumberOfA, _ := strconv.Atoi(statistic.NumberOfA)
		NumberOfC, _ := strconv.Atoi(statistic.NumberOfC)
		NumberOfG, _ := strconv.Atoi(statistic.NumberOfG)
		NumberOfT, _ := strconv.Atoi(statistic.NumberOfT)
		AverageQuality, _ := strconv.ParseFloat(statistic.AverageQuality, 32)
		averageQualityTotal = ((float64(averageQualityTotal) * float64(numberOfNucleotidesTotal)) + (float64(AverageQuality) * (float64(NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)))) / (float64(numberOfNucleotidesTotal) + float64(NumberOfN+NumberOfA+NumberOfC+NumberOfG+NumberOfT))
		numberOfNucleotidesTotal = numberOfNucleotidesTotal + (NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)
		if _, ok := numberOfReadsMap[statistic.Channel]; ok {
			numberOfReadsMap[statistic.Channel] += NumberOfReads
		} else {
			numberOfReadsMap[statistic.Channel] = NumberOfReads
		}
		if _, ok := averageQualityMap[statistic.Channel]; ok {
			averageQualityMap[statistic.Channel] = ((averageQualityMap[statistic.Channel] * float64(numberOfNucleotidesMap[statistic.Channel])) + (float64((NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)) * AverageQuality)) / (float64(numberOfNucleotidesMap[statistic.Channel] + NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT))
		} else {
			averageQualityMap[statistic.Channel] = AverageQuality
		}
		if _, ok := numberOfNucleotidesMap[statistic.Channel]; ok {
			numberOfNucleotidesMap[statistic.Channel] += (NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)
		} else {
			numberOfNucleotidesMap[statistic.Channel] = (NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)
		}
		if _, ok := numberOfNMap[statistic.Channel]; ok {
			numberOfNMap[statistic.Channel] += NumberOfN
		} else {
			numberOfNMap[statistic.Channel] = NumberOfN
		}
		if _, ok := numberOfAMap[statistic.Channel]; ok {
			numberOfAMap[statistic.Channel] += NumberOfA
		} else {
			numberOfAMap[statistic.Channel] = NumberOfA
		}
		if _, ok := numberOfCMap[statistic.Channel]; ok {
			numberOfCMap[statistic.Channel] += NumberOfC
		} else {
			numberOfCMap[statistic.Channel] = NumberOfC
		}
		if _, ok := numberOfGMap[statistic.Channel]; ok {
			numberOfGMap[statistic.Channel] += NumberOfG
		} else {
			numberOfGMap[statistic.Channel] = NumberOfG
		}
		if _, ok := numberOfTMap[statistic.Channel]; ok {
			numberOfTMap[statistic.Channel] += NumberOfT
		} else {
			numberOfTMap[statistic.Channel] = NumberOfT
		}
		if _, ok := maxReadMap[statistic.Channel]; ok {
			if (NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT) > maxReadMap[statistic.Channel] {
				maxReadMap[statistic.Channel] = (NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)
			}
		} else {
			maxReadMap[statistic.Channel] = (NumberOfN + NumberOfA + NumberOfC + NumberOfG + NumberOfT)
		}*/
		if _, ok := rawDataLengthMap[statistic.ChannelNumber]; ok {
			rawDataLengthMap[statistic.ChannelNumber] += RawDataLength
		} else {
			rawDataLengthMap[statistic.ChannelNumber] = RawDataLength
		}
	}
}

func init() {
	metrics.UseNilMetrics = true
	prometheus.MustRegister(version.NewCollector("fast5_exporter"))
}

func toFlagString(name string, help string, value string) *string {
	flag.CommandLine.String(name, value, help) // hack around flag.Parse and klog.init flags
	return kingpin.Flag(name, help).Default(value).String()
}

func toFlagBool(name string, help string, value bool, valueString string) *bool {
	flag.CommandLine.Bool(name, value, help) // hack around flag.Parse and klog.init flags
	return kingpin.Flag(name, help).Default(valueString).Bool()
}

func toFlagStringsVar(name string, help string, value string, target *[]string) {
	flag.CommandLine.String(name, value, help) // hack around flag.Parse and klog.init flags
	kingpin.Flag(name, help).Default(value).StringsVar(target)
}

func toFlagStringVar(name string, help string, value string, target *string) {
	flag.CommandLine.String(name, value, help) // hack around flag.Parse and klog.init flags
	kingpin.Flag(name, help).Default(value).StringVar(target)
}

func toFlagBoolVar(name string, help string, value bool, valueString string, target *bool) {
	flag.CommandLine.Bool(name, value, help) // hack around flag.Parse and klog.init flags
	kingpin.Flag(name, help).Default(valueString).BoolVar(target)
}

func toFlagIntVar(name string, help string, value int, valueString string, target *int) {
	flag.CommandLine.Int(name, value, help) // hack around flag.Parse and klog.init flags
	kingpin.Flag(name, help).Default(valueString).IntVar(target)
}

func main() {
	var (
		ontFast5DirPath = toFlagString("ont-fast5-dir-path", "Path to the dir where fast5 from ONT sequencer will be stored.", "/tmp")
		listenAddress   = toFlagString("web.listen-address", "Address to listen on for web interface and telemetry.", ":9307")
		metricsPath     = toFlagString("web.telemetry-path", "Path under which to expose metrics.", "/metrics")
		topicFilter     = toFlagString("topic.filter", "Regex that determines which topics to collect.", ".*")
		groupFilter     = toFlagString("group.filter", "Regex that determines which consumer groups to collect.", ".*")
		logSarama       = toFlagBool("log.enable-sarama", "Turn on Sarama logging, default is false.", false, "false")

		opts = exporterOpts{}
	)

	toFlagStringsVar("fast5.server", "Address (host:port) of fast5 server.", "fast5:9092", &opts.uri)
	toFlagBoolVar("sasl.enabled", "Connect using SASL/PLAIN, default is false.", false, "false", &opts.useSASL)
	toFlagBoolVar("sasl.handshake", "Only set this to false if using a non-fast5 SASL proxy, default is true.", true, "true", &opts.useSASLHandshake)
	toFlagStringVar("sasl.username", "SASL user name.", "", &opts.saslUsername)
	toFlagStringVar("sasl.password", "SASL user password.", "", &opts.saslPassword)
	toFlagStringVar("sasl.mechanism", "The SASL SCRAM SHA algorithm sha256 or sha512 or gssapi as mechanism", "", &opts.saslMechanism)
	toFlagStringVar("sasl.service-name", "Service name when using kerberos Auth", "", &opts.serviceName)
	toFlagStringVar("sasl.kerberos-config-path", "Kerberos config path", "", &opts.kerberosConfigPath)
	toFlagStringVar("sasl.realm", "Kerberos realm", "", &opts.realm)
	toFlagStringVar("sasl.kerberos-auth-type", "Kerberos auth type. Either 'keytabAuth' or 'userAuth'", "", &opts.kerberosAuthType)
	toFlagStringVar("sasl.keytab-path", "Kerberos keytab file path", "", &opts.keyTabPath)
	toFlagBoolVar("sasl.disable-PA-FX-FAST", "Configure the Kerberos client to not use PA_FX_FAST, default is false.", false, "false", &opts.saslDisablePAFXFast)
	toFlagBoolVar("tls.enabled", "Connect to fast5 using TLS, default is false.", false, "false", &opts.useTLS)
	toFlagStringVar("tls.server-name", "Used to verify the hostname on the returned certificates unless tls.insecure-skip-tls-verify is given. The fast5 server's name should be given.", "", &opts.tlsServerName)
	toFlagStringVar("tls.ca-file", "The optional certificate authority file for fast5 TLS client authentication.", "", &opts.tlsCAFile)
	toFlagStringVar("tls.cert-file", "The optional certificate file for fast5 client authentication.", "", &opts.tlsCertFile)
	toFlagStringVar("tls.key-file", "The optional key file for fast5 client authentication.", "", &opts.tlsKeyFile)
	toFlagBoolVar("server.tls.enabled", "Enable TLS for web server, default is false.", false, "false", &opts.serverUseTLS)
	toFlagBoolVar("server.tls.mutual-auth-enabled", "Enable TLS client mutual authentication, default is false.", false, "false", &opts.serverMutualAuthEnabled)
	toFlagStringVar("server.tls.ca-file", "The certificate authority file for the web server.", "", &opts.serverTlsCAFile)
	toFlagStringVar("server.tls.cert-file", "The certificate file for the web server.", "", &opts.serverTlsCertFile)
	toFlagStringVar("server.tls.key-file", "The key file for the web server.", "", &opts.serverTlsKeyFile)
	toFlagBoolVar("tls.insecure-skip-tls-verify", "If true, the server's certificate will not be checked for validity. This will make your HTTPS connections insecure. Default is false", false, "false", &opts.tlsInsecureSkipTLSVerify)
	toFlagStringVar("fast5.version", "Fast5 version", sarama.V2_0_0_0.String(), &opts.fast5Version)
	toFlagBoolVar("use.consumelag.zookeeper", "if you need to use a group from zookeeper, default is false", false, "false", &opts.useZooKeeperLag)
	toFlagStringsVar("zookeeper.server", "Address (hosts) of zookeeper server.", "localhost:2181", &opts.uriZookeeper)
	toFlagStringVar("fast5.labels", "fast5 cluster name", "", &opts.labels)
	toFlagStringVar("refresh.metadata", "Metadata refresh interval", "30s", &opts.metadataRefreshInterval)
	toFlagBoolVar("offset.show-all", "Whether show the offset/lag for all consumer group, otherwise, only show connected consumer groups, default is true", true, "true", &opts.offsetShowAll)
	toFlagBoolVar("concurrent.enable", "If true, all scrapes will trigger fast5 operations otherwise, they will share results. WARN: This should be disabled on large clusters. Default is false", false, "false", &opts.allowConcurrent)
	toFlagIntVar("topic.workers", "Number of topic workers", 100, "100", &opts.topicWorkers)
	toFlagBoolVar("fast5.allow-auto-topic-creation", "If true, the broker may auto-create topics that we requested which do not already exist, default is false.", false, "false", &opts.allowAutoTopicCreation)
	toFlagIntVar("verbosity", "Verbosity log level", 0, "0", &opts.verbosityLogLevel)

	plConfig := plog.Config{}
	plogflag.AddFlags(kingpin.CommandLine, &plConfig)
	kingpin.Version(version.Print("fast5_exporter"))
	kingpin.HelpFlag.Short('h')
	kingpin.Parse()

	labels := make(map[string]string)

	if opts.labels != "" {
		for _, label := range strings.Split(opts.labels, ",") {
			splitted := strings.Split(label, "=")
			if len(splitted) >= 2 {
				labels[splitted[0]] = splitted[1]
			}
		}
	}

	setup(*ontFast5DirPath, *listenAddress, *metricsPath, *topicFilter, *groupFilter, *logSarama, opts, labels)
}

func setup(
	ontFast5DirPath string,
	listenAddress string,
	metricsPath string,
	topicFilter string,
	groupFilter string,
	logSarama bool,
	opts exporterOpts,
	labels map[string]string,
) {
	klog.InitFlags(flag.CommandLine)
	if err := flag.Set("logtostderr", "true"); err != nil {
		klog.Errorf("Error on setting logtostderr to true: %v", err)
	}
	err := flag.Set("v", strconv.Itoa(opts.verbosityLogLevel))
	if err != nil {
		klog.Errorf("Error on setting v to %v: %v", strconv.Itoa(opts.verbosityLogLevel), err)
	}
	defer klog.Flush()

	klog.V(INFO).Infoln("Starting fast5_exporter", version.Info())
	klog.V(DEBUG).Infoln("Build context", version.BuildContext())

	totalSizeMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "total", "size"),
		"Stats",
		[]string{"address", "name"}, labels,
	)
	totalAverageQualityMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "total", "average_quality"),
		"Stats",
		[]string{"address", "name"}, labels,
	)
	numberOfReadsMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "amount", "reads"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	nMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "amount", "N"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	aMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "amount", "A"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	cMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "amount", "C"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	gMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "amount", "G"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	tMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "amount", "T"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	maxReadMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "max", "read"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	averageQualityMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "average", "quality"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)
	rawDataLengthMetric = prometheus.NewDesc(
		prometheus.BuildFQName(namespace, "raw_data", "length"),
		"Stats",
		[]string{"address", "name", "channel"}, labels,
	)

	if logSarama {
		sarama.Logger = log.New(os.Stdout, "[sarama] ", log.LstdFlags)
	}

	exporter, err := NewExporter(opts, topicFilter, groupFilter)
	if err != nil {
		//			klog.Fatalln(err)
	}
	//		defer exporter.client.Close()
	prometheus.MustRegister(exporter)

	http.Handle(metricsPath, promhttp.Handler())
	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		_, err := w.Write([]byte(`<html>
		        <head><title>Fast5 Exporter</title></head>
		        <body>
		        <h1>Fast5 Exporter</h1>
		        <p><a href='` + metricsPath + `'>Metrics</a></p>
		        </body>
		        </html>`))
		if err != nil {
			klog.Error("Error handle / request", err)
		}
	})
	http.HandleFunc("/healthz", func(w http.ResponseWriter, r *http.Request) {
		_, err := w.Write([]byte("ok"))
		if err != nil {
			klog.Error("Error handle /healthz request", err)
		}
	})

	if opts.serverUseTLS {
		klog.V(INFO).Infoln("Listening on HTTPS", listenAddress)

		_, err := CanReadCertAndKey(opts.serverTlsCertFile, opts.serverTlsKeyFile)
		if err != nil {
			klog.Error("error reading server cert and key")
		}

		clientAuthType := tls.NoClientCert
		if opts.serverMutualAuthEnabled {
			clientAuthType = tls.RequireAndVerifyClientCert
		}

		certPool := x509.NewCertPool()
		if opts.serverTlsCAFile != "" {
			if caCert, err := ioutil.ReadFile(opts.serverTlsCAFile); err == nil {
				certPool.AppendCertsFromPEM(caCert)
			} else {
				klog.Error("error reading server ca")
			}
		}

		tlsConfig := &tls.Config{
			ClientCAs:                certPool,
			ClientAuth:               clientAuthType,
			MinVersion:               tls.VersionTLS12,
			CurvePreferences:         []tls.CurveID{tls.CurveP521, tls.CurveP384, tls.CurveP256},
			PreferServerCipherSuites: true,
			CipherSuites: []uint16{
				tls.TLS_ECDHE_RSA_WITH_AES_256_GCM_SHA384,
				tls.TLS_ECDHE_RSA_WITH_AES_128_GCM_SHA256,
				tls.TLS_ECDHE_RSA_WITH_AES_256_CBC_SHA,
				tls.TLS_ECDHE_RSA_WITH_AES_128_CBC_SHA256,
				tls.TLS_RSA_WITH_AES_256_GCM_SHA384,
				tls.TLS_RSA_WITH_AES_256_CBC_SHA,
				tls.TLS_RSA_WITH_AES_128_CBC_SHA256,
			},
		}
		server := &http.Server{
			Addr:      listenAddress,
			TLSConfig: tlsConfig,
		}
		klog.Fatal(server.ListenAndServeTLS(opts.serverTlsCertFile, opts.serverTlsKeyFile))
	} else {
		klog.V(INFO).Infoln("Listening on HTTP", listenAddress)
		klog.Fatal(http.ListenAndServe(listenAddress, nil))
	}
}
