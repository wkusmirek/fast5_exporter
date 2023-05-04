import argparse
import json
import os
from collections import Counter
import numpy as np

def get_file_size(file_path):
    """
    Returns the size of the filepath in bytes.
    :param filePath: The path to file in the watch directory
    :return:
    """
    return os.path.getsize(file_path)

def check_is_pass(path, avg_quality):
    """
    Check that the fastq file is contains pass or fail reads
    Parameters
    ----------
    path
    avg_quality

    Returns
    -------

    """
    folders = os.path.split(path)
    if "pass" in folders[0]:
        return True
    elif "fail" in folders[0]:
        return False
    else:
        if avg_quality != None:
            if avg_quality >= 7:
                return True  # This assumes we have been unable to find either pass
                # or fail and thus we assume the run is a pass run.
            else:
                return False
        else:
            return True

def average_quality(quality):
    """
    Return the average quality of a read from it's quality string
    Parameters
    ----------
    quality: str

    Returns
    -------
    float
    """
    return -10 * np.log10(
        np.mean(
            10
            ** (
                np.array(
                    (-(np.array(list(quality.encode("utf8"))) - 33)) / 10, dtype=float
                )
            )
        )
    )

def parse_fastq_description(description):
    """
    Parse the description found in a fastq reads header

    Parameters
    ----------
    description: str
        A string of the fastq reads description header

    Returns
    -------
    description_dict: dict
        A dictionary containing the keys and values found in the fastq read headers
    """
    description_dict = {}
    descriptors = description.split(" ")
    # Delete symbol for header
    for item in descriptors:
        if "=" in item:
            bits = item.split("=")
            description_dict[bits[0]] = bits[1]
    return description_dict

def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield desc, name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield desc, name, seq, "".join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield desc, name, seq, None  # yield a fasta record instead
                break

def split_all(path):
    """
    Split the path into it's relative parts so we can check all the way down the tree
    :param path: Path provided by the user to the watch directory
    :return: A list of paths created from the base path, so we can use relative and absolute paths
    """
    all_parts = []
    while 1:
        parts = os.path.split(path)

        if parts[0] == path:  # sentinel for absolute paths
            all_parts.insert(0, parts[0])
            break
        elif parts[1] == path:  # sentinel for relative paths
            all_parts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            all_parts.insert(0, parts[1])
    return all_parts

def check_fastq_path(path):
    """
    Check if the files we are looking at are directly inside pass or fail folders, or one or two removed.
    :param path:
    :return:
    """
    folders_in_path = split_all(path)
    try:
        if folders_in_path[-2] in ("pass", "fail"):
            return "{}_{}".format(folders_in_path[-2], folders_in_path[-1])
        elif folders_in_path[-3] in ("pass", "fail"):
            return "{}_{}_{}".format(folders_in_path[-3], folders_in_path[-2], folders_in_path[-1])
        else:
            return "{}".format(folders_in_path[-1])
    except Exception as e:
        return "{}".format(folders_in_path[-1])

def get_runid(fastq):
    """
    Open a fastq file, read the first line and parse out the Run ID
    :param fastq: path to the fastq file to be parsed
    :type fastq: str
    :return runid: The run ID of this fastq file as a string
    """
    runid = ""
    if fastq.endswith(".gz"):
        with gzip.open(fastq, "rt") as file:
            for _ in range(1):
                line = file.readline()
    else:
        with open(fastq, "r") as file:
            for _ in range(1):
                line = file.readline()
    for _ in line.split():
        if _.startswith("runid"):
            runid = _.split("=")[1]
    return runid

def parse_fastq_file(
    fastq_path #, run_dict, args, header, minotour_api, sequencing_stats
):
    """
    
    Parameters
    ----------
    fastq_path: str
        Path to the fastq file to parse
    run_dict: dict
        Dictionary containing the run
    args: argparse.NameSpace
        The namespace of the parsed command line args
    header: dict
        The authentication header to make requests to the server
    minotour_api: minFQ.minotourapi.MinotourAPI
        The minotour API class instance for convenience request to the server
    sequencing_stats: minFQ.utils.SequencingStatistics
        The class to track the sequencing upload metrics and files

    Returns
    -------
    int 
        The updated number of lines we have already seen from the unblocked read ids file
    """
    numberOfReads = 0
    numberOfN = 0
    numberOfA = 0
    numberOfC = 0
    numberOfG = 0
    numberOfT = 0
    run_dict = dict()
    # The updated number of lines we have already seen from the unblocked read ids file
    #log.debug("Parsing fastq file {}".format(fastq_path))
    # Get runId from the path
    run_id = get_runid(fastq_path)
    #print(run_id)
    payload = {
        "file_name": str(check_fastq_path(fastq_path)),
        "run_id": run_id,
        "md5": "0",
        "run": None
    }
    #if run_id in run_dict:
    #    fastq_file = minotour_api.post(
    #        EndPoint.FASTQ_FILE, json=payload, base_id=run_id
    #    )
    #if not sequencing_stats.fastq_info[run_id]["run_id"]:
    #    sequencing_stats.fastq_info[run_id]["run_id"] = run_id
    counter = 0
    quality = ''
    handle = gzip.open if fastq_path.endswith(".gz") else open
    with handle(fastq_path, "rt") as fh:
        for desc, name, seq, qual in readfq(fh):
            description_dict = parse_fastq_description(desc)
            counter += 1
            numberOfReads += 1
            numberOfN = seq.count('N') + seq.count('n')
            numberOfA = seq.count('A') + seq.count('a')
            numberOfC = seq.count('C') + seq.count('c')
            numberOfG = seq.count('G') + seq.count('g')
            numberOfT = seq.count('T') + seq.count('t')
            quality += qual
            #sequencing_stats.reads_seen += 1
            #sequencing_stats.fastq_info[run_id]["reads_seen"] += 1
            #sequencing_stats.fastq_message = "processing read {}".format(counter)
            try:
                ### We haven't seen this run before - so we need to check stuff.
                # create runCOllection in the runDict
                if description_dict["runid"] not in run_dict:
                    create_run_collection(
                        run_id,
                        run_dict,
                        args,
                        header,
                        description_dict,
                        sequencing_stats,
                    )
                    fastq_file = minotour_api.post(
                        EndPoint.FASTQ_FILE, json=payload, base_id=run_id
                    )
                    run_dict[run_id].get_readnames_by_run(fastq_file["id"])
                parse_fastq_record(
                    name,
                    seq,
                    qual,
                    fastq_path,
                    run_dict,
                    args,
                    fastq_file,
                    counter=counter,
                    sequencing_statistic=sequencing_stats,
                    description_dict=description_dict,
                )

            except Exception as e:
                #sequencing_stats.reads_corrupt += 1
                continue
    # This chunk of code will mean we commit reads every time we get a new file?
    #for runs in run_dict:
    #    run_dict[runs].read_names = []
    #    run_dict[runs].commit_reads()
    #try:
        # Update the fastq file entry on the server that says we provides file size to database record
    #    payload = {
    #        "file_name": str(check_fastq_path(fastq_path)),
    #        "run_id": run_id,
    #        "md5": get_file_size(fastq_path),
    #        "run": run_dict[run_id].run,
    #    }
    #    fastq_file = minotour_api.put(EndPoint.FASTQ_FILE, json=payload, base_id=run_id)
    #except Exception as err:
    #    log.error("Problem with uploading file {}".format(err))
    #sequencing_stats.time_per_file = time.time()
    #sequencing_stats.fastq_info[run_id]["files_processed"] += 1
    averageQuality = round(average_quality(quality), 2,)
    isPass = check_is_pass(fastq_path, averageQuality)
    fileSize = get_file_size(fastq_path)
    read = description_dict.get("read", None)
    runId = description_dict.get("runid", None)
    channel = description_dict.get("ch", None)
    startTime = description_dict.get("start_time", None)
    return numberOfReads, numberOfN, numberOfA, numberOfC, numberOfG, numberOfT, averageQuality, isPass, fileSize, read, runId, channel, startTime

def main():
    parser = argparse.ArgumentParser(
        description="Parse fastq file"
    )
    parser.add_argument(
        "--path", default="", help="Path to fastq file"
    )
    args = parser.parse_args()

    numberOfReads, numberOfN, numberOfA, numberOfC, numberOfG, numberOfT, averageQuality, isPass, fileSize, read, runId, channel, startTime = parse_fastq_file(args.path)

    result = list()
    vars = dict()
    vars['numberOfReads'] = str(numberOfReads)
    vars['numberOfN'] = str(numberOfN)
    vars['numberOfA'] = str(numberOfA)
    vars['numberOfC'] = str(numberOfC)
    vars['numberOfG'] = str(numberOfG)
    vars['numberOfT'] = str(numberOfT)
    vars['averageQuality'] = str(averageQuality)
    vars['isPass'] = str(isPass)
    vars['fileSize'] = str(fileSize)
    vars['read'] = str(read)
    vars['runId'] = str(runId)
    vars['channel'] = str(channel)
    vars['startTime'] = str(startTime)
    result.append(vars)

    print('[' + ','.join(json.dumps(d) for d in result) + ']')

if __name__ == "__main__":
    main()
