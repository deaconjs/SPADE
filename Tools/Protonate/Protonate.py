
def __get_hydrogens (pchain,chain_name, res_number, atom_name):
    """ Return a list of hydrogens that are associated with an atom
    """
    hydAtoms=[]
    hydrogen_pat=re.compile("^d*H")
    target_res = pchain.residues_dict[str(res_number)]
    # Loop through all of the atoms and search for the hydrogens
    for atom in target_res.atoms:
        # If it is a hydrogen atom then find its closest non hydrogen
        # atom and that atom is its owner. If the owner is atom_name
        # then add it to the list of hydAtoms
        if (hydrogen_pat.search(atom.atom_type)):
            hyd=atom
            minDist=-1
            owner=None
            for atom2 in target_res.atoms:
                if (not hydrogen_pat.search(atom2.atom_type) and
                    (hyd.dist(atom2) < minDist or minDist == -1)):
                    owner = atom2
                    minDist = hyd.dist(atom2)
            # If the closest atom is the atom_name then add it to the hydrogen
            # list
            if (owner.atom_type == target_res.atoms_dict[atom_name].atom_type):
                hydAtoms.append (hyd)
    return hydAtoms

def __find_atoms_for_protonation (pchain,protonsInfo):
    """ Find all of the atoms that can be protonated in this protein.
        Returns a list of protons where an element contains
        {'atom','aa'[,'prev_aa'],'protonInfo'}.
    """
    protons = []
    for j in range(len(pchain.residues)):
        aa = pchain.residues[j]
    
        # Add all of the amino acid specific protons
        for i in range(len(aa.atoms)):
            # Find the proton information for this atom if it has any
            protonInfo = find_proton_info (aa.res_type,aa.atoms[i],protonsInfo)
            if (protonInfo != None):
                protons.append ({'atom': aa.atoms[i], 'aa': aa, 'protonInfo': protonInfo})

        # if it is the first amino acid in the chain (N-TERMINUS)
        if (j == 0):
            protonInfo = find_proton_info ('N-TERMINUS',aa.atoms_dict['N'],protonsInfo)
            if (protonInfo != None):
                protons.append ({'atom': aa.atoms_dict['N'], 'aa': aa, 'protonInfo': protonInfo})
        # Add the backbone amino acid, which is common to all amino acids except
        # the first amino acid in a chain
        elif (aa.res_type != 'PRO'): # Every other residue except PRO
            protonInfo = find_proton_info ('BACKBONE',aa.atoms_dict['N'],protonsInfo)
            if (protonInfo != None):
                # Store the previous amino acid because we need its carbon
                protons.append ({'atom': aa.atoms_dict['N'], 'aa': aa, 'prev_aa': pchain.residues[j-1],
                                 'protonInfo': protonInfo})
    return protons

def protonate (pchain,protonFile,redo=False):
    """ Protonate a protein. Add the hydrogens to the protein. i.e. create a list
        of possible donors and all of their hydrogens 'atomsHydList': {'atom','hydAtoms'}
    """
    pchain.atomsHydList=[]

    protonsInfo = read_protonation_info (protonFile)
    protons = pchain.__find_atoms_for_protonation (protonsInfo)

    # Position the hydrogens
    for proton in protons:
        D=proton['atom']
        aa=proton['aa']
        hydAtoms=[]
        # Initialize the hydrogen atoms to any that are already in the protein
        hydAtoms = pchain.__get_hydrogens(D.chain_name,D.res_number,D.atom_type)
        if (len(hydAtoms) != 0 and redo == False):
            # Save the list of hyrdogen atoms for this atom 
            pchain.atomsHydList.append ({'atom': D, 'hydAtoms': hydAtoms})
        elif (len(hydAtoms) != 0 and redo == True):
            print ("TODO: Add deletion of the current hydrogens")
  
        # sp2, 1H 2DD
        if (proton['protonInfo']['hyb'] == 'sp2' and
            proton['protonInfo']['bonds'] == '1H 2DD'):
            DD1Name=proton['protonInfo']['DD1Name']
            angle_offset=0 # 0 -> DD1-D-H = DD2-D-H
            if (DD1Name == 'None'):
                print "ERROR: No DD1 atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DD1=aa.atoms_dict[DD1Name]
            if (proton['protonInfo']['angles'] == 'DD1-D-H = DD2-D-H'):
                DD2Name=proton['protonInfo']['DD2Name']
                if (DD2Name=='None'):
                    print "ERROR: No DD2 atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD2=aa.atoms_dict[DD2Name]
            elif (proton['protonInfo']['angles'] == '(C-N-H)-(CA-N-H)=4; C CA N H are planar'):
                prev_aa=proton['prev_aa']
                DD2Name=proton['protonInfo']['DD2Name']
                if (DD2Name=='None'):
                    print "ERROR: No DD2 atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD2=prev_aa.atoms_dict[DD2Name]
                angle_offset = 4
            # Find the angle C-N-H to place the hydrogen
            Dv=r_[D.x,D.y,D.z]
            # Make Dv the origin
            DD2v=r_[DD2.x,DD2.y,DD2.z]-Dv
            DD1v=r_[DD1.x,DD1.y,DD1.z]-Dv
            Dv=r_[0,0,0]
            theta=acos(dot(DD1v,DD2v)/(mag(DD1v)*mag(DD2v)))
            angle=pi-theta/2+radians(angle_offset)
            hydPos=findPlanarPosition (proton['protonInfo']['D-H'],angle*180/pi,D,DD2,DD1,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos))
      
        # sp2, 1H 1DD
        elif (proton['protonInfo']['hyb'] == 'sp2' and
              proton['protonInfo']['bonds'] == '1H 1DD'):
            DDName=proton['protonInfo']['DD1Name']
            if (DDName == 'None'):
                print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DD=aa.atoms_dict[DDName]
            DDDName=proton['protonInfo']['DDD1Name']
            if (DDDName=='None'):
                print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DDD=aa.atoms_dict[DDDName]
            # This configuration has two mutually exclusive hydrogen positions
            # TODO: For now we are placing both hydrogens but we need to add some kind of
            # collision resolution
            # TODO: Change find planar to take an array of angles...
            hydPos1=findPlanarPosition (proton['protonInfo']['D-H'],250,D,DD,DDD,True)
            hydPos2=findPlanarPosition (proton['protonInfo']['D-H'],110,D,DD,DDD,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1))
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos2))

        # sp2, 2H 1DD
        elif (proton['protonInfo']['hyb'] == 'sp2' and
              proton['protonInfo']['bonds'] == '2H 1DD'):
            DDName=proton['protonInfo']['DD1Name']
            if (DDName == 'None'):
                print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DD=aa.atoms_dict[DDName]
            DDDName=proton['protonInfo']['DDD1Name']
            if (DDDName=='None'):
                print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DDD=aa.atoms_dict[DDDName]
            # This configuration has two hydrogen positions
            # TODO: Change find planar to take an array of angles...
            hydPos1=findPlanarPosition (proton['protonInfo']['D-H'],120,D,DD,DDD,True)
            hydPos2=findPlanarPosition (proton['protonInfo']['D-H'],360-120,D,DD,DDD,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1)) 
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos2))

        # sp3, 1H 1DD
        elif (proton['protonInfo']['hyb'] == 'sp3' and
              proton['protonInfo']['bonds'] == '1H 1DD'):
            DDName=proton['protonInfo']['DD1Name']
            if (DDName == 'None'):
                print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DD=aa.atoms_dict[DDName]
            DDDName=proton['protonInfo']['DDD1Name']
            if (DDDName=='None'):
                print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DDD=aa.atoms_dict[DDDName]
            if (proton['protonInfo']['angles'] == 'DD-D-H=110'):
                hydPos1=findCirclePosition (proton['protonInfo']['D-H'],110,0,D,DD,DDD,True)
            elif (proton['protonInfo']['angles'] == 'DD-D-H=96'):
                hydPos1=findCirclePosition (proton['protonInfo']['D-H'],96,0,D,DD,DDD,True)
            elif (proton['protonInfo']['angles'] == 'DD-D-H=110; Circle Angle=240'):
                hydPos1=findCirclePosition (proton['protonInfo']['D-H'],110,240,D,DD,DDD,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1))

        # sp3, 3H 1DD
        elif (proton['protonInfo']['hyb'] == 'sp3' and
              proton['protonInfo']['bonds'] == '3H 1DD'):
            DDName=proton['protonInfo']['DD1Name']
            if (DDName == 'None'):
                print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DD=aa.atoms_dict[DDName]
            DDDName=proton['protonInfo']['DDD1Name']
            if (DDDName=='None'):
                print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                continue
            DDD=aa.atoms_dict[DDDName]
            hydPos1=findCirclePosition (proton['protonInfo']['D-H'],110,0,D,DD,DDD,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1))
            hydPos2=findCirclePosition (proton['protonInfo']['D-H'],110,120,D,DD,DDD,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos2))
            hydPos3=findCirclePosition (proton['protonInfo']['D-H'],110,240,D,DD,DDD,True)
            hydAtoms.append(pchain.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos3))

        # Save the list of hyrdogen atoms for this atom 
        pchain.atomsHydList.append ({'atom': D, 'hydAtoms': hydAtoms})
    
    pchain.parent.renumber_atoms()

def get_avail_donors (pchain,donorFile):
    donorsInfo=read_donor_info (donorFile)

    # Find all of the atoms that have hydrogens and that can donate
    availDonors=[]
    for atomHyd in pchain.atomsHydList:
        donorInfo = find_donor_info (atomHyd['atom'].parent.res_type,atomHyd['atom'],donorsInfo)
        if (donorInfo != None):
            availDonors.append ({'donorAtom': atomHyd['atom'], 'hydAtoms': atomHyd['hydAtoms']})

    return availDonors

def get_avail_acceptors (pchain,accFile):
    accsInfo=read_acc_info(accFile)

     # Find all of the atoms in the protein that can be acceptors
    availAcceptors = []
    for j in range(len(pchain.residues)):
        aa = pchain.residues[j]
        for i in range(len(aa.atoms)):
            # Find the acceptor information if it has any
            accInfo = find_acc_info (aa.res_type,aa.atoms[i],accsInfo)
            if (accInfo != None):
                availAcceptors.append ({'accAtom': aa.atoms[i],#'accInfo': accInfo, Don't need this?
                                        'AAAtom': aa.atoms_dict[accInfo['AAName']]})
    return availAcceptors

def find_as_neighborhood (pchain,as_residue_inxs):
    """ Find the neighborhood of amino acids that has the largest number
        of active site residues
        Return: a list of residues
    """
    maxInx=-1
    maxCnt=-1
    # Find the neighborhood with the maximum number of as residues
    for i in range(len(as_residue_inxs)):
        as_rex = as_residue_inxs[i]
        cnt=0
        for rex in range(len(pchain.residues[as_rex].neighbors)):
            if (rex in as_residue_inxs): # if this neighbor is an as residue
                cnt+=1
        if (maxCnt < cnt):
            maxInx=i
            maxCnt=cnt

    if (maxInx < 0):
        print "ERROR: incorrect index (find_as)"

    # Create the list of residues based on the neighborhood selected
    as_residues = [pchain.residues[as_residue_inxs[maxInx]]]
    for rex in pchain.residues[as_residue_inxs[maxInx]].neighbors:
        if (rex in as_residue_inxs):
            as_residues.append(pchain.residues[rex])
    return as_residues

def find_as_concat (pchain,as_residue_inxs):
    """ Find the active site of the protein by grouping residues based on whether
        a residue can see at least one of the other residues in the group
        Return: a list of residues
    """
    graph=[]
    as_residues=[]
    # Create the graph representation of the neighborhoods
    for i in range(len(as_residue_inxs)):
        rex_row=as_residue_inxs[i]
        row = []
        for j in range(len(as_residue_inxs)):
            rex_col=as_residue_inxs[j]
            if (j == i): # Of course we can see ourselves
                row.append (1)
            # Can see each other
            elif (rex_col in pchain.residues[rex_row].neighbors):
                row.append(1)
            else: # Cannot see each other
                row.append(0)
        graph.append(row)

    """
    # Test case
    # Should have three groups:
    # [0, 4, 3, 6, 1]
    # [2, 5, 9]
    # [7, 8]
    graph=[[1,0,0,0,1,0,0,0,0,0],
           [0,1,0,0,0,0,1,0,0,0],
           [0,0,1,0,0,1,0,0,0,0],
           [0,0,0,1,1,0,0,0,0,0],
           [1,0,0,1,1,0,1,0,0,0],
           [0,0,1,0,0,1,0,0,0,1],
           [0,1,0,0,1,0,1,0,0,0],
           [0,0,0,0,0,0,0,1,1,0],
           [0,0,0,0,0,0,0,1,1,0],
           [0,0,0,0,0,1,0,0,0,1]]
    """

    # Show the graph
    print 'Initial Graph'
    for row in graph:
        print row

    # Keep track of which residues have been attached to
    # a group
    attached=[]
    for i in range(len(graph)):
        attached.append(0)
        
    # Perform a BFT to get the possible active sites
    groups = []
    # While there is still a residue that is not attached
    while (0 in attached):
        # Find a residue that has not been attached
        for z in range(len(attached)):
            if (attached[z]==0):
                break
        group = [z]
        attached[z]=1
        stack=[]
        stack.append(z)

        # Traverse the graph
        while (len(stack) > 0):
            row = graph[stack.pop()]
            for j in range(len(row)):
                # Connected the current node and not already attached
                # to a group
                if (row[j] == 1 and attached[j]==0): 
                    group.append(j)
                    attached[j]=1
                    stack.append(j)
        groups.append (group)

    # Print the resultant groups
    print 'Final Groups'
    for group in groups:
        print group
        
    return as_residues


def find_as_concat_n (pchain,n,as_residue_inxs):
    """ Find the active site of the protein by grouping residues based on whether
        a residue can see at least 'n' of the other residues in the group
        Return: a list of residues
    """
    graph=[]
    as_residues=[]
    # Create the graph representation of the neighborhoods
    for i in range(len(as_residue_inxs)):
        rex_row=as_residue_inxs[i]
        row = []
        for j in range(len(as_residue_inxs)):
            rex_col=as_residue_inxs[j]
            if (j == i): # Of course we can see ourselves
                row.append (1)
            # Can see each other
            elif (rex_col in pchain.residues[rex_row].neighbors):
                row.append(1)
            else: # Cannot see each other
                row.append(0)
        graph.append(row)

    """
    # Test case
    # Should have three groups:
    # [0, 4, 3, 6, 1]
    # [2, 5, 9]
    # [7, 8]
    graph=[[1,0,0,0,1,0,0,0,0,0],
           [0,1,0,0,0,0,1,0,0,0],
           [0,0,1,0,0,1,0,0,0,0],
           [0,0,0,1,1,0,0,0,0,0],
           [1,0,0,1,1,0,1,0,0,0],
           [0,0,1,0,0,1,0,0,0,1],
           [0,1,0,0,1,0,1,0,0,0],
           [0,0,0,0,0,0,0,1,1,0],
           [0,0,0,0,0,0,0,1,1,0],
           [0,0,0,0,0,1,0,0,0,1]]
    """

    # Show the graph
    print 'Initial Graph'
    for row in graph:
        print row

    # Keep track of which residues have been attached to
    # a group
    attached=[]
    for i in range(len(graph)):
        attached.append(0)
        
    # Perform a BFT to get the possible active sites
    groups = []
    # While there is still a residue that is not attached
    while (0 in attached):
        # Find a residue that has not been attached
        for head in range(len(attached)):
            if (attached[head]==0):
                break
        group = [head]
        attached[head]=1
        stack=[]
        stack.append(head)

        # Traverse the graph
        while (len(stack) > 0):
            row = graph[stack.pop()]
            for j in range(len(row)):
                # Connected to the current node and not already attached
                # to a group
                if (row[j] == 1 and attached[j]==0): 
                    group.append(j)
                    attached[j]=1
                    stack.append(j)
        groups.append (group)

    # Find out if each amino acid in the group can see at least two
    # of the other amino acids by examining the neighborhoods
    for i in range(len(groups)):
        group=groups[i]
        # For each member of the group count the number it can
        # see
        group_tmp=[]
        for inx in group:
            residue=pchain.residues[as_residue_inxs[inx]]
            cnt=0
            # Search through all of the other members of the group
            for inx2 in group:
                if (inx == inx2): continue
                # If it is the neighborhood then increment the count
                if (as_residue_inxs[inx2] in residue.neighbors):
                    cnt+=1
            if (cnt >= n):
                group_tmp.append(inx)
        print 'Before:',group
        print 'After: ',group_tmp
        groups[i]=group_tmp
                

    # Print the resultant groups
    print 'Final Groups'
    for group in groups:
        print group

    # Find the group with the largest number of residues = active site
    maxInx=-1
    maxCnt=-1
    for i in range(len(groups)):
        group = groups[i]
        if (len(group) > maxCnt):
            maxCnt=len(group)
            maxInx=i

    # Using that group get the actual residue numbers
    if (maxInx < 0):
        print "ERROR: incorrect index (find_as)"

    # Create the list of residues based on the group selected
    as_residues = []
    for inx in groups[maxInx]:
        as_residues.append(pchain.residues[as_residue_inxs[inx]])

    return as_residues
