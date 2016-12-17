from ftplib import FTP
import string
import Tkinter
import Pmw
import os
import gzip
import shutil

def _pdb_callback(list, item):
	if item != 'total':
		list.append(item)
		
def _execute_rcsb_download(ftp, list, window, target_directory):
	# make sure the directory exists
	if not os.path.isdir(target_directory):
		os.mkdir(target_directory)
	# get the items to download
	list = window.pdb_combo.getvalue()
	for request in list:
		final_name = request[3:7]
		# I happen to have some malicious code that changes extentions of downloaded filenames
		# and heres the fix
		if request[8:11] != 'ent':
			request = request[:8] + 'ent' + request[11:]
		final_path = os.path.join(target_directory, final_name)
		# make sure there's a directory for the pdb file
		if not os.path.isdir(final_path):
			os.mkdir(final_path)
		final_file = os.path.join(final_path, final_name + '.pdb.Z')
		# retrieve the file 
		ftp.retrbinary('RETR ' + request, open(final_file, 'wb').write)
	window.destroy()
	ftp.quit()
	
def download_from_rcsb(target_directory):
	""" target_directory should be a directory (present or not) within /Databases. Include
	only the path from there. """
	# transform target_directory into its full path
	target_directory = os.path.join(os.path.abspath('.'), '.', target_directory)
	# acquire the source directory
	ftp = FTP('ftp.rcsb.org', 'anonymous', '')   # connect to host, default port
	ftp.cwd("pub/pdb/data/structures/all/pdb")
	available_pdbs = []
	c_lambda = lambda item: _pdb_callback(available_pdbs, item)
	ftp.retrlines('NLST', c_lambda)
	# make a request window with a combo box for selection
	tnsf_win = Tkinter.Toplevel()
	tnsf_win.pdb_combo = Pmw.ScrolledListBox(tnsf_win,
											 label_text = 'Select files for download:',
											 labelpos = 'n',
											 items = available_pdbs,
											 listbox_selectmode='multiple')
	tnsf_win.pdb_combo.pack(fill = 'x', expand = 1, padx = 8, pady = 8)
	c_lambda = lambda: _execute_rcsb_download(ftp, list, tnsf_win, target_directory)
	tnsf_win.ok_button = Tkinter.Button(tnsf_win, text='OK', command=c_lambda)
	tnsf_win.ok_button.pack()

def _execute_database_transfer(tnsf_win, target_directory):
	# Create a set of directories in the target directory
	# Move over the files, renaming .ent files
	if not os.path.isdir(target_directory):
		os.mkdir(target_directory)
	source_path = os.path.join(os.path.abspath('.'), tnsf_win.source_entry.getvalue())
	files = os.listdir(source_path)
	for file in files:
		move_source = os.path.join(source_path, file)
		if file[-3:] == 'ent':
			if not os.path.isdir(os.path.join(target_directory, file[3:7])):
				os.mkdir(os.path.join(target_directory, file[3:7]))
			move_target = os.path.join(target_directory, file[3:7], file[3:7] + '.pdb')
		elif file[-3:] == 'pdb' and len(string.split(file,'_')) <= 2:
			if not os.path.isdir(os.path.join(target_directory, file[:-4])):
				os.mkdir(os.path.join(target_directory, file[:-4]))
			move_target = os.path.join(target_directory, file[:-4], file[:-4] + '.pdb')
		shutil.move(move_source, move_target)
	
def transform_from_pdbs(target_directory):
	# Build a form to acquire a source directory
	tnsf_win = Tkinter.Toplevel()
	# create an absolute path for the target
	target_directory = os.path.join(os.path.abspath('.'), 'Databases', target_directory)
	# Add an entry for the directory and an ok button
	tnsf_win.source_entry = Pmw.EntryField(tnsf_win, labelpos = 'w', label_text = 'Target Database (from /SPADE)', validate = None, value='default_database')
	tnsf_win.source_entry.pack(side='top',fill='x', expand=1, padx=10, pady=5)
	c_lambda = lambda: _execute_database_transfer(tnsf_win, target_directory)
	tnsf_win.ok_button = Tkinter.Button(tnsf_win, text='OK', command=c_lambda)
	tnsf_win.ok_button.pack(side='top')
	


