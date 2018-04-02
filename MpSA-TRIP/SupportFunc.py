#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
import os, csv, regex, datetime, json, logging, logging.config
from toolshed import nopen
"""
**** SUPPORT FUNCTIONS
"""

def GetTotalSeqRecords(input_file):
    '''
    This function count of number strings in fastq file and return sequences number
    '''
    with nopen(input_file) as f:
        TotalSeqRecords = int(sum(1 for _ in f)) / 4
    return TotalSeqRecords

def head(list_or_dict, start=0, n=10):
    start, n = int(start), int(int(start) + int(n))
    if type(list_or_dict) == dict:
        for i in dict(list_or_dict.items()[start:n]): print i, dict(list_or_dict.items()[start:n])[i]
    elif type(list_or_dict) == list:
        for i in list_or_dict[start:n]: print i

def simpleWrite(list_or_dict, path, name, delim="\t"):
    filePath = os.path.join(path, name)
    with open(filePath, "wb") as handle:
        if type(list_or_dict) == dict:
            header = ["id", "value"]
            writer = csv.DictWriter(handle, fieldnames=header, delimiter=delim)
            writer.writeheader()
            for i in list_or_dict:
                writer.writerow({
                    "id": i,
                    "value": list_or_dict[i]})
        elif type(list_or_dict) == list:
            header = ["id"]
            writer = csv.DictWriter(handle, fieldnames=header, delimiter=delim)
            writer.writeheader()
            for i in list_or_dict:
                writer.writerow({
                    "id": i})

def jsonWrite(var_dict, outpath, outfilename):
    path_to_out = os.path.join(outpath, outfilename)
    with open(path_to_out, 'wb') as handle:
        json.dump(var_dict, handle)

def EstimateCalculationTime(obj):
    if type(obj) == dict or type(obj) == list:
        l = float(len(obj))
        velocity = 5*10**5
        sec = (float(l**2)/2)/velocity
        now = datetime.datetime.now().strftime("%H:%M")
        if sec <= 60:
            t = str(float(sec)) + " sec. Current time: " + now
        elif sec <= 3600 and sec > 60:
            t = str(float(sec)/(60)) + " min. Current time: " + now
        elif sec <= 24*60**2 and sec > 3600:
            t = str(float(sec)/(60**2)) + "h."
        else:
            t = str(float(sec)/(24*60**2)) + "d. Current time: " + now
        return t
    return "I don't estimate calculation time for this object. Sorry..."

def makeStatFromBowtieAlign(bwtAlignerDict):
    LogInfo('Some statistics from bowtie aligner')
    expr = regex.compile('[ ]*(?P<count>\d.*) ((\((?P<pct>[0-9\.\%].*)\) .*aligned (?P<times>.*time.*|))|reads.*)')
    for key in bwtAlignerDict:
        LogInfo("Bowtie align report for {}".format(key))
        for line in bwtAlignerDict[key]:
            m = expr.match(line)
            if m is not None:
                if m.group('times') is not None:
                    LogInfo("     Aligned {}: {} ({})".format(m.group('times'), m.group('count'), m.group('pct')))
                else:
                    LogInfo("     {} input reads".format(m.group('count')))


def setup_logging(default_path='logging.json', default_level=logging.INFO, env_key='LOG_CFG'):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)

def LogInfo(msg):
    logger = logging.getLogger(__name__)
    logger.info(msg)

def LogErr(msg):
    logger = logging.getLogger(__name__)
    logger.error(msg)