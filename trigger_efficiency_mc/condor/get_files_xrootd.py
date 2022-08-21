import subprocess, os, re

primary_name='CRAB_PrivateMC_RunII_UL_2017'
dataset = 'Jpsi_50to100_Dstar_DPS_13TeV' #Jpsi_50to100_Dstar_DPS_13TeV
crab_folder = '220422_173510' #220422_173510
n_folders = 9

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


cmds = [ 
   f'xrdfs xrootd-redir.ultralight.org ls -u /store/group/uerj/mabarros/{primary_name}/{dataset}/{crab_folder}/{i:04}/' for i in range(n_folders)
]
#print(cmds)

cat = ""
out_file = dataset + "_path.txt"
for i in cmds:
   os.system(f"{i} > {i.split('/')[-2]}")
   cat += f" {i.split('/')[-2]}"

os.system(f"cat {cat} > {out_file}")
os.system(f"rm -rf {cat}")
file_list = open(out_file, 'r').readlines()

for idx, f in enumerate(file_list):
   file_list[idx] = re.sub('transfer-\d*', 'xrootd-redir', f)

list(set(file_list))
file_list.sort(key=natural_keys)

final_list = []
for i in file_list:
    if i not in final_list:
        final_list.append(i)

with open(out_file, 'w') as f:
    for i in final_list:
        f.write(i)
