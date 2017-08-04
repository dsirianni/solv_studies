import sys

class geom():

    def __init__(self, xyz, solvent = ''):
        #Assumes xyz is list of lists or xyz file
        if type(xyz) == str: 
            geometry = []
            with open(xyz,'r') as xyz_fil:
                next(xyz_fil)
                next(xyz_fil)
                for lin in xyz_fil:
                    if len(lin.split()) != 4:
                        print("Ill formatted xyz.")
                        sys.exit()
                    name = lin.split()[0]
                    coords = [float(x) for x in lin.split()[1:]]
                    geometry.append([name,coords[0],coords[1],coords[2]])
            self.fil = xyz.split('.xyz')[0]
        else: 
            geometry = xyz
            self.fil = "MOLECULE"
        self.geometry = geometry
        name = ''
        for atom in geometry: name+=atom[0].upper()
        self.name = name
        solvent = solvent.upper()
        self.solvent = solvent
    
    def same_name(self, other_name):
        lother_name = list(other_name)
        lname = list(self.name)
        for l in self.name: 
            try: 
                lother_name.remove(l)
                lname.remove(l)
            except: pass
        if len(lother_name) == 0 and len(lname) == 0: same = True
        else: same = False
        return same 

    def dist(self,a,b):
        xdist = (a[1]-b[1])**2
        ydist = (a[2]-b[2])**2
        zdist = (a[3]-b[3])**2
        dist = (xdist+ydist+zdist)**0.5
        return dist
    
    def intersect(self,frag,frags):
        intersects = set() 
        for f in frags:
            if frag & f == set([]): continue
            intersects = intersects | frag | f
            frags.remove(f) 
        return intersects
    
    def bound(self,a,b):
        cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
          'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
          'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
          'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
          'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
          'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
          'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}
        if cov_rad[a[0]]+cov_rad[b[0]] >=  self.dist(a,b): return True
        else: return False
    
    def bfs(self):
        minifrags = []
        #make sets of bonded atoms
        for i,si in enumerate(self.geometry):
            for j,sj in enumerate(self.geometry):
                if j>=i: continue
                if self.bound(si,sj): minifrags.append(set([str(si),str(sj)]))
        
        frags = []
        while minifrags != []:
            one = self.intersect(minifrags[0],minifrags)
            last = []
            while one != set([]):
                one = self.intersect(one,minifrags)
                if one != set([]): last = list(one)
            if last != []:
                listlast = []
                for atom in last:
                    #messy conversion from set of strings to list of atoms 
                    name = atom.split("'")[1]
                    coords = atom.split(',')[1:]
                    coords[-1] = coords[-1].split(']')[0]
                    coords = [float(x) for x in coords]
                    listlast.append([name,coords[0],coords[1],coords[2]])
                frags.append(geom(listlast))
        self.frags = frags
        return frags

    def write_input(self, file_name = '', solv_mon = 'C'):
        self.monsets()
        solv_mon = solv_mon.upper()
        if file_name == '': file_name = self.fil + '.in'
        with open(file_name,'w') as inp_fil:            
            inp_fil.write("molecule "+self.fil+" {\n")
            inp_fil.write("0 1\n")
            for a in self.frags[0]: inp_fil.write(a[0]+' '+' '.join([str(x) for x in a[1:]])+'\n') 
            if len(self.frags[2]) != 0 and solv_mon == 'A':
                for c in self.frags[2]: 
                    for catom in c: inp_fil.write(catom[0]+' '+' '.join([str(x) for x in catom[1:]])+'\n') 
            inp_fil.write("--\n0 1\n")
            for b in self.frags[1]: inp_fil.write(b[0]+' '+' '.join([str(x) for x in b[1:]])+'\n') 
            if len(self.frags[2]) != 0 and solv_mon == 'B':
                for c in self.frags[2]: 
                    for catom in c: inp_fil.write(catom[0]+' '+' '.join([str(x) for x in catom[1:]])+'\n') 
            if len(self.frags[2]) != 0 and solv_mon == 'C':
                inp_fil.write("--\n0 1\n")
                for c in self.frags[2]: 
                    for catom in c: inp_fil.write(catom[0]+' '+' '.join([str(x) for x in catom[1:]])+'\n') 
            inp_fil.write('''units angstrom\nno_com\nno_reorient\nsymmetry c1
}\n\nmemory 30000 MB\n
set {
  basis         jun-cc-pvdz
  df_basis_scf  jun-cc-pvdz-jkfit
  df_basis_sapt jun-cc-pvdz-ri
  freeze_core   true
  scf_type df
  guess sad
  minao_basis cc-pvtz-minao
  maxiter 100
  ints_tolerance 0.0
}\n\nenergy('fisapt0')\n\n''')


    def monsets(self):
        self.bfs()
        names = []
        #assumes there will only be three sets:
        #solvent molecules, monA, monB
        fragA = []
        fragB = []
        fragC = []
        for frag in self.frags:
            if frag.same_name(self.solvent): fragC.append(frag.geometry)
            elif fragA == []: fragA = frag.geometry
            elif fragB == []: fragB = frag.geometry
            else: 
                print("Make sure there are only three unique molecules.")
                sys.exit()
        self.frags = [fragA,fragB,fragC]




