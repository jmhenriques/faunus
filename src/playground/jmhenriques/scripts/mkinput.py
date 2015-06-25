#!/usr/bin/env python
#
import sys, math, operator, json, ast
#
"""
Masses are correct for the 'free' aminoacids.
Protonated forms should have +1 in mass, but that is irrelevant.

All pKa's taken from Evers_2012_Langmuir_28_11843.
Except the pKa for O-phospho-L-serine (pSE), taken from Zachariou_1996_JPC_100_12680.
"""
#
__author__  = "Joao Henriques"
__email__   = "joao.henriques@teokem.lu.se"
__date__    = "2015.06.25"
__status__  = "Production"
#
"""
Key:
resn (str) --- charge (int) --- mass (float) --- hydrophobic (str) --- pKa (float)
"""
dict={'ALA' :[ 0, 71 , "True" ],
      'ARG' :[ 0, 156, "False", 12.5],
      'HARG':[ 1, 156, "False"],
      'ASN' :[ 0, 114, "False"],
      'ASP' :[-1, 115, "False", 3.9 ],
      'HASP':[ 0, 115, "False"],
      'CYS' :[-1, 103, "False", 8.5 ],
      'HCYS':[ 0, 103, "False"],
      'GLU' :[-1, 129, "False", 4.1 ],
      'HGLU':[ 0, 129, "False"],
      'GLN' :[ 0, 128, "False"],
      'GLY' :[ 0, 57 , "False"],
      'HIS' :[ 0, 137, "False", 6.5 ],
      'HHIS':[ 1, 137, "False"],
      'ILE' :[ 0, 113, "True" ],
      'LEU' :[ 0, 113, "True" ],
      'LYS' :[ 0, 128, "False", 10.8],
      'HLYS':[ 1, 128, "False"],
      'MET' :[ 0, 131, "True" ],
      'nSE' :[ 0, 165, "False"],
      'PHE' :[ 0, 147, "True" ],
      'PRO' :[ 0, 97 , "True" ],
      'pSE' :[-2, 165, "False", 5.8 ],
      'HpSE':[-1, 165, "False"],
      'SER' :[ 0, 87 , "False"],
      'THR' :[ 0, 101, "False"],
      'TRP' :[ 0, 186, "True" ],
      'TYR' :[-1, 163, "False", 10.1],
      'HTYR':[ 0, 163, "False"],
      'VAL' :[ 0, 99 , "True" ],
      'CTR' :[-1, 16 , "False", 3.6 ],
      'HCTR':[ 0, 16 , "False"],
      'NTR' :[ 0, 14 , "False", 8.6 ],
      'HNTR':[ 1, 14 , "False"]}
#
def Usage():
    print "\nTo produce a JSON file: \n\t%s json <input fasta file> <pH>" %(sys.argv[0])
    print "\nTo produce a AAM file : \n\t%s aam  <input fasta file> <bond length (Angstrom)>\n" %(sys.argv[0])
#
def readFasta(file):
    """
    Reads a FASTA sequence and returns a list with the corresponding 
    three-letter code.
    Appends NTR and CTR to both ends. X and Z don't exist, but I added 
    them for pSE and nSE.
    """
    list = ['NTR']
    dict = {'A' : 'ALA',
            'R' : 'ARG',
            'N' : 'ASN',
            'D' : 'ASP',
            'C' : 'CYS',
            'E' : 'GLU',
            'Q' : 'GLN',
            'G' : 'GLY',
            'H' : 'HIS',
            'I' : 'ILE',
            'L' : 'LEU',
            'K' : 'LYS',
            'M' : 'MET',
            'F' : 'PHE',
            'P' : 'PRO',
            'S' : 'SER',
            'T' : 'THR',
            'W' : 'TRP',
            'Y' : 'TYR',
            'V' : 'VAL',
            'X' : 'pSE',
            'Z' : 'nSE'
            }
    f = open(file, 'r')
    for line in f.readlines():
        if line[0] != '>':
            for elem in line:
                if elem.upper() in dict:
                    list.append(dict[elem])
    list.append('CTR')
    return list
#
def getDensity(fasta):
    """
    Calculates the protein density using an equation 
    taken from Fischer_2004_ProteinScience_13_2825.
    """
    mass = float(0)
    for res in fasta:
        if res in dict:
            mass += dict[res][1]/1000.0           # kDa
    density = 1.41 + 0.145 * math.exp(-mass/13.0) # g/cm^3
    return density
#
def getRadius(mass, density):
    """
    Computes the radius of an aa residue with a given (molecular) volume,
    assuming a spherical approximation. The (molecular) volume is obtained 
    from the residue mass and protein density.
    """
    NA = 6.0221413e+23                      # mol^-1
    molarVol = float(mass)/float(density)   # cm^3/mol
    volume = (molarVol*math.pow(1e8, 3))/NA # Angstrom^3
    r = math.pow(((3.0*volume)/(4.0*math.pi)), 1.0/3.0)
    return r
#
def writeJSON(pH, density):
    """
    Simplest way of writing a fully automatic JSON file.
    """
    vals1 = {}
    vals2 = {}
    for key in dict:
        val = dict[key]
        if len(key)<4 and len(val)==4:
            vals1.update({"H-"+key:{"bound":"H"+key, "free":key, "pKd":val[3], "pX":round(float(pH), 1)}})
        vals2.update({key:{"q":val[0], "mw":round(val[1], 1), "r":round(getRadius(val[1], density), 1), "hydrophobic":ast.literal_eval(val[2])}})
    data = {"processes":vals1, "atomlist":vals2}
    print json.dumps(data, indent=4, sort_keys=True)
#            
def writeAAM(fasta, density, bondlen):
    """
    Writes an AAM file given a three letter a.a. sequence.
    I need to write a FASTA reader for Faunus, because this is a dumb
    way of dealing with this.
    """
    print len(fasta)
    c=1
    i=0.0
    for res in fasta:
        if res in dict:
            # AAM key: 
            # resn x y z charge mass radius
            print "%s\t%3i %6.1f 0.0 0.0 %2i %5.1f %5.1f" %(res, c, i, dict[res][0], dict[res][1], getRadius(dict[res][1], density))
        c += 1
        i += float(bondlen)
#
if __name__=="__main__":
    if len(sys.argv[:]) > 1:
        if sys.argv[1].lower() in ("h", "-h", "--h", "help", "-help", "--help"):
            Usage()
        elif sys.argv[1].lower()=='json':
            fname = sys.argv[2]
            fasta = readFasta(fname)
            density = getDensity(fasta)
            pH = sys.argv[3]
            writeJSON(pH, density)
        elif sys.argv[1].lower()=="aam":
            fname = sys.argv[2]
            fasta = readFasta(fname)
            density = getDensity(fasta)
            bondlen = sys.argv[3]
            writeAAM(fasta, density, bondlen)
        else:
            Usage()
    else:
        Usage()
