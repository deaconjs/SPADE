
from Bio import SCOP
from os import *
import parms
import ScopDomainViewer

class ScopIterator:
	#isAprofile: 1 (if scopStructure is from a pre-saved profile) or 0 (if it is the original scopStructure)
	def __init__(self, isAprofile):
		self.ScopStructure = None
		self.ScopDomains = None
		self.currentDomain = -1
		self.currentDomainLocation = ""
		
		self.isAprofile = isAprofile
		self.profileName = ""  #will only be set if the scop structure is a profile
		self.profileSunid = "" #same as above

	#Load the scop structure
	def loadScop(self, cla_path=None, des_path=None, hie_path=None):
		if cla_path is None:
			self.cla_file = file(parms.get('scopClassification'), 'r')
			self.des_file = file(parms.get('scopDescription'), 'r')
			self.hie_file = file(parms.get('scopHierarchy'), 'r')
		elif cla_path is not None and hie_path is not None and des_path is not None:
			self.cla_file = file(cla_path, 'r')
			self.des_file = file(des_path, 'r')
			self.hie_file = file(hie_path, 'r')
		else:
			print "Error in arguments"
			return

		self.ScopStructure = SCOP.Scop(self.cla_file, self.des_file, self.hie_file)
		self.ScopDomains = self.ScopStructure.getDomains()
		
		cla_file_name = self.cla_file.name
		if self.isAprofile == 1:
			self.profileName = cla_file_name[cla_file_name[:cla_file_name.rindex('\\')].rindex('\\')+1:cla_file_name.rindex('\\')]
			self.profileSunid = cla_file_name[cla_file_name.rindex('\\')+1:cla_file_name.rindex('\\')+6]
			print "Profile Name: "+self.profileName
			
					
	def load_first_domain(self):
		if self.ScopStructure is None or self.ScopDomains is None:
			print "SCOP Structure not loaded"
			return
		
		self.currentDomain = 0
		self.currentDomainLocation = self.get_current_domain_location()
	
	def load_next_domain(self):
		if self.currentDomain == -1:
			self.load_first_domain()
		else:
			self.currentDomain = self.currentDomain + 1
			self.currentDomainLocation = self.get_current_domain_location()
		
	#Returns the location of the domain and if does not exist, returns ""
	#Location determined based on self.isAprofile.
	def get_current_domain_location(self):	
		if self.currentDomain == -1:
			print "Load the first domain first before accessing this funciton"
		
		if self.isAprofile == 1:
			profilePath = parms.get('profilesPath') + self.profileName
			tmpPath = ""
			#parse the domain's parent structure to get the path		
			curDomain = self.currentDomainNode()
			hierarchy = []
			while curDomain is not None:
				name = self.filterFoldername(curDomain.description)
				hierarchy.append(name)
				if self.profileSunid.find(curDomain.sunid) == 0:
					curDomain = None	
				else:
					curDomain = curDomain.parent 	
			
			while len(hierarchy) > 0:
				foldr = hierarchy.pop()
				profilePath = profilePath +'\\'+foldr

			if not path.isfile(profilePath+'\\'+self.currentDomainNode().sid+'.ent'):
				print "Cannot find File!"
				return ""
			return profilePath+'\\'+self.currentDomainNode().sid+'.ent'
		
		else:
			pdb_path = parms.get('scopPDBsLocation') + self.currentDomainNode().sid[2:4] + '\\' + self.currentDomainNode().sid + '.ent'
			if path.isfile(pdb_path):
				return pdb_path
	
	#returns lineage as a list
	def getCurrentDomainLineage(self):
		type = parms.get('classifications_fullname')
		lineage = []

		currentNode = self.currentDomainNode()
		while currentNode is not None and len(currentNode.type) > 0:
			lineage.insert(0, type[currentNode.type]+":  "+currentNode.description+'\n')
			currentNode = currentNode.parent
	
		if self.isAprofile == 1:
			self.hie_file.seek(0)
			line = self.hie_file.readline()
			while line[0] == '#':
				lineage.insert(0,line[1:].strip()+'\n')
				line = self.hie_file.readline()
				
		return lineage
	
	def currentDomainNode(self):
		return self.ScopDomains[self.currentDomain]
	
	#Removes illegal characters from directory names
	def filterFoldername(self, name):
		dirNamesDonotIncludelist = parms.get('dirNamesDonotIncludelist')	
		for itm in dirNamesDonotIncludelist:
			name = name.replace(itm, '')
		return name.strip()	
	
	def getDomainBySid(self, sid):
		pos = 0
		for item in self.ScopDomains:
			if item.sid.find(sid) == 0:
				self.currentDomain = pos
				break
			pos = pos + 1

	def getDomainBySunid(self, sunid):
		pos = 0
		for item in self.ScopDomains:
			if item.sunid.find(sunid) == 0:
				self.currentDomain = pos
				break
			pos = pos + 1	
			
	def getDomainsListByPdbid(self, pdbid):
		pos = 0
		domsList = []
		for item in self.ScopDomains:
			if item.residues.pdbid.find(pdbid) == 0:
				domsList.append(item)
			pos = pos + 1				
		return domsList
	
	#The following "are****" functions take in a list of domains and return true if
	#they belong to one type(class, family, etc)
	def areDomainsOfSameClass(self, doms):
		type = doms[0].sccs[0]
		for item in doms:
			if item.sccs[0] != type:
				return 0
		return 1
	def areDomainsOfSameFold(self, doms):
		type = doms[0].sccs[:3]
		for item in doms:
			if item.sccs.find(type) != 0:
				return 0
		return 1		
	def areDomainsOfSameSuperFamily(self, doms):
		type = doms[0].sccs[:5]
		for item in doms:
			if item.sccs.find(type) != 0:
				return 0
		return 1		
		
	def areDomainsOfSameFamily(self, doms):	
		type = doms[0].sccs[:7]
		for item in doms:
			if item.sccs.find(type) != 0:
				return 0
		return 1				
		
if __name__ == '__main__':

	prof = ScopIterator(1)
	cla_file = "C:\\CourseWork\\CS499\\MyTest\\SCOPSearchProfiles\\1101010\\46457_cla.txt"
	des_file = "C:\\CourseWork\\CS499\\MyTest\\SCOPSearchProfiles\\1101010\\46457_des.txt"
	hie_file = "C:\\CourseWork\\CS499\\MyTest\\SCOPSearchProfiles\\1101010\\46457_hie.txt"

	prof.loadScop(cla_file, des_file, hie_file)
	prof.load_first_domain()
	prof.load_next_domain()
	print prof.getCurrentDomainLineage()
	prof.getDomainBySunid('85678')
	prof.getDomainBySid('d1ngkf_')
	prof.getDomainsListByPdbid('1ngk')
	prof.areDomainsOfSameFamily(prof.getDomainsListByPdbid('1ngk'))
	
	