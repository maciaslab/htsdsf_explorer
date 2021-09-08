import glob,os

nohide =["webapp","data","plateinfo","persistent","README.md","HTSDSF.exe","settings.ini"]
for file in glob.glob("*"):
    if file in nohide:
        print ("NOHIDE",file)
    else:
        print("HIDE",file)
        proc=os.popen('attrib +h ' + file)
        proc.read()
        proc.close()