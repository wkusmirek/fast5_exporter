import schedule
import time
import shutil
import os

files = os.listdir('data')

def job():
    shutil.copyfile('data/'+files[0], '/tmp/fast5/'+files[0])
    files.pop(0)

schedule.every(10).seconds.do(job)

while True:
    schedule.run_pending()
    time.sleep(1)
