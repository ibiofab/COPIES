import psutil, time, io, zipfile, json, os, sys, base64, glob, traceback, contextlib, shutil, random, subprocess
os.environ["LAMBDA_LIBS"] = "/mnt/efs/venv38/lib/python3.8/site-packages/"
os.environ["SEARCH_IDX"] = "/mnt/efs/data_new/search_idx.txt"

#sys.path.append(os.environ["LAMBDA_LIBS"])
import boto3
dynamodb = boto3.resource('dynamodb', region_name='us-east-2')
dynamo_table = dynamodb.Table('ibiofab-copies')
    
def read_param(p_name):
    with open("/sys/fs/cgroup/{}".format(p_name), "r") as f:
        return int(f.read())

def event_handler(event):
        

    working_directory = "/tmp/output_%010x" % random.randrange(16**10)
    os.makedirs(working_directory)   
    
    command = ["python", "/mnt/efs/COPIES/code/main.py"]
    for key in event["body"].keys():
        # TODO: Sanitize input.
        command.append("--"+str(key))
        command.append(str(event["body"][key]))

    command.extend(["-out", working_directory+"/output_data.csv"])

    with open(os.path.join(working_directory,"output.txt"), "w") as f:
        #pass
        start = time.perf_counter()
        sub = subprocess.Popen(command, stdout=f, stderr=f)
        sub.wait()
        sub2 = subprocess.Popen(["python", "/mnt/efs/COPIES/code/create_viz.py", os.path.join(working_directory,"output_data.csv"),
                                os.path.join(working_directory,"copies_visualization.html")])
        sub2.wait()
        mem_peak = int(read_param("memory/memory.max_usage_in_bytes")/(1024*1024))
        mem_lim = int(read_param("memory/memory.limit_in_bytes")/(1024*1024))
        vcpu_share = read_param("cpu/cpu.cfs_quota_us")/read_param("cpu/cpu.cfs_period_us")
        p = subprocess.run("curl http://169.254.169.254/latest/meta-data/instance-type", shell=True, capture_output=True)
        instance_type = p.stdout.decode("utf-8")
        f.write("=================================\nCommand Used\n")
        f.write(" ".join(command))
        f.write("\n=================================")
        f.write("\nvCPUs, Mem Limit (MiB), Execution time, Mem Peak (Mib), instance_type\n")
        f.write("{},{},{},{},{}".format(psutil.cpu_count(), mem_lim, time.perf_counter()-start, mem_peak, instance_type))
    

    memfile = io.BytesIO()
    memzip = zipfile.ZipFile(memfile, mode="x", compression=zipfile.ZIP_DEFLATED, compresslevel=9)

    old_dir = os.getcwd()
    os.chdir(working_directory)
    for fname in [x for x in glob.glob("**", recursive=True) if os.path.isfile(x)]:
        memzip.write(fname)
    memzip.close()
    os.chdir(old_dir)
    memfile.seek(0)
    
    s3 = boto3.client('s3', region_name='us-east-2' )#, config=botocore.config.Config(s3={'addressing_style':'path'}))
    s3.upload_fileobj(memfile, "ibiofab-copies", event["entry_id"], ExtraArgs={"ContentDisposition": "attachment; filename=copies_Result.zip"})
    try:
        s3.upload_file(os.path.join(working_directory,"copies_visualization.html"), "ibiofab-copies", event["entry_id"]+"_graph" , ExtraArgs={"ContentType" : "text/html", "ContentDisposition": "inline"})
    except FileNotFoundError:
        # generation script did not run correctly.
        pass

    dynamo_table.update_item(Key={"id":str(event["entry_id"])},
    UpdateExpression="SET percent_complete = :q",
    ExpressionAttributeValues={":q" : 100} )

        
    return

event_handler(json.loads(sys.argv[1]))
