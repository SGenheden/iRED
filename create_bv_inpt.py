# Author: Samuel Genheden samuel.genheden@gmail.com, 2011-2015

"""
Program to create a cpptraj input for iRED analysis
"""

import argparse

import def_bv

if __name__ == '__main__':

  disclaimer="""
If -p and -t is not supplied on the command-line, the program will ask for them
\n
If a question mark (?) is specified as the vector tag, all available vector tags
 will be printed to the screen.
  """
  # Parse the command-line
  parser = argparse.ArgumentParser(description="Program to create cpptraj input for iRED analysis",epilog=disclaimer,formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-v', '--vectors',help="the vectors to analyse",default="?")
  parser.add_argument('-p', '--pdbfile',help="a template PDB file")
  parser.add_argument('-t', '--trajectory',nargs='+',help="a list of trajectories to analyse")
  parser.add_argument('-o', '--output',help="the name of the cpptraj output",default="iredout")
  args = parser.parse_args()

  if args.vectors == "?" :
    print "The available bond vector types are: %s"%(" ".join(def_bv.BondVectors.vectorkeys()))
    quit()

  # Ask for the essential input options if they were not supplied on the command-line
  if args.pdbfile == None :
    print "Enter the name of a template PDB file: ",
    args.pdbfile = raw_input().strip()
  if args.trajectory == None :
    print "Enter the name of trajectory file to analyse: ",
    args.trajectory = [raw_input().strip()]   

  # Create and initialise a class that will parse the bond vectors
  vectors = def_bv.BondVectors.vectors(args.vectors)(args.pdbfile)
  
  # Produce a cpptraj input file and dump it to standard output
  print "parm %s"%args.pdbfile
  for s in args.trajectory :
    print "trajin ",s
  for i,vec in enumerate(vectors) :
    print "vector v%d %s"%(i+1,vec)
  print "matrix ired name matired order 2"
  print "diagmatrix matired vecs %d name ired.vec"%len(vectors)
  print "ired order 2 modes ired.vec orderparamfile %s out %s"%(args.output+"_s2",args.output)
