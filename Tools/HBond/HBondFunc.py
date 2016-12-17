import re
import string

"""
 Read the protonation information from the data file
"""
def read_protonation_info (protonFile):
    protonsInfo = []
    # Read the protonation information from the data file
    fproton = open (protonFile,"r")
    comment_pat = re.compile ("^#")
    # Format of a protonation line is:
    # AA Name,D,DD1,DD2,DDD1,DDD2,Hybridization,Bonds,D-H distance,Angles
    for line in fproton:
        line = string.strip (line)
        if (not comment_pat.search (line) and line != ''):
            fields=string.split (line,",")
            protonInfo={}
            protonInfo['aminoAcidName']=fields[0]
            protonInfo['atomName']=fields[1]
            protonInfo['DD1Name']=fields[2]
            protonInfo['DD2Name']=fields[3]
            protonInfo['DDD1Name']=fields[4]
            protonInfo['DDD2Name']=fields[5]
            protonInfo['hyb']=fields[6]
            protonInfo['bonds']=fields[7]
            protonInfo['D-H']=string.atof(fields[8])
            protonInfo['angles']=fields[9]
            protonsInfo.append (protonInfo)
    fproton.close()

    return protonsInfo

"""
 Find the proton information for a given atom and residue
"""
def find_proton_info (res_type,atom,protonsInfo):
    for protonInfo in protonsInfo:
        atom_type = atom.atom_type
        if (protonInfo['atomName'] == atom_type and
            (res_type == protonInfo['aminoAcidName'] or protonInfo['aminoAcidName'] == '*')):
            return protonInfo
    return None

"""
 Find the donor information for a given atom
"""
def find_donor_info (res_type,atom,donorsInfo):
    for donorInfo in donorsInfo:
        if (donorInfo['atomName'] == atom.atom_type and
            (res_type == donorInfo['aminoAcidName'] or donorInfo['aminoAcidName'] == '*')):
            return donorInfo
    return None

"""
 Read donor information
"""
def read_donor_info (donorFile):
    donorsInfo = []
    fdonor = open (donorFile,"r")
    comment_pat = re.compile ("^#")
    # Format of a donor line is:
    # AA Donor Name, Donor Atom Name
    for line in fdonor:
        line = string.strip (line)
        if (not comment_pat.search (line) and line != ''):
            fields=string.split (line,",")
            donorInfo={}
            donorInfo['aminoAcidName']=fields[0]
            donorInfo['atomName']=fields[1]
            donorsInfo.append (donorInfo)
    fdonor.close()
    return donorsInfo


"""
 Find the acceptor information for a given atom
"""
def find_acc_info (res_type,atom,accsInfo):
    for accInfo in accsInfo:
        if (accInfo['atomName'] == atom.atom_type and
            (res_type == accInfo['aminoAcidName'] or accInfo['aminoAcidName'] == '*')):
            return accInfo
    return None

"""
 Read the acceptor information
"""
def read_acc_info (accFile):
    accsInfo = []
    facc = open (accFile,"r")
    comment_pat = re.compile ("^#")
    # Format of a acceptor line is:
    # Amino Acid Name, A, AA
    for line in facc:
        line = string.strip (line)
        if (not comment_pat.search (line) and line != ''):
            fields=string.split (line,",")
            accInfo={}
            accInfo['aminoAcidName']=fields[0]
            accInfo['atomName']=fields[1]
            accInfo['AAName']=fields[2]
            accsInfo.append (accInfo)
    facc.close()
    return accsInfo
