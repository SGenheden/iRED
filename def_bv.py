# Author: Samuel Genheden samuel.genheden@gmail.com, 2011-2015

"""
Module containing classes for parsing iRED bond vectors from a PDB file
"""

import inspect
import sys

class BondVector(object) :
  """
  Class to store state of single bond vector
  
  Attributes
  ----------
  resname : string
    the residue name
  resid : string
    the residue serial number
  atom1 : string
    the name of the first atom
  atom2 : string
    the name of the second atom
  """
  def __init__(self,resnam,resid,atom1,atom2) :
    self.resnam = resnam
    self.resid = resid
    self.atom1 = atom1
    self.atom2 = atom2
  def __str__(self) :
    return ":%s@%s ired :%s@%s"%(self.resid,self.atom1,self.resid,self.atom2)

class BondVectors(object) :
  """
  Base class for classes that parse PDB files to define iRED bond vectors
  
  Each subclass should implement the "_parse routine" to correctly parse the PDB file
  and define a class variable "tag" to specify a short string used to select the class.
  
  The BondVectors class is an iterator class and iterating over it returns the defined
  bond vectors. 
  
  The len() operator is implemented and returns the number of defined vectors.
  
  It defines two class methods:
  vectorkeys : list of string
    the tag of the defined subclasses
  vectors : class
    the subclass of BondVectors, selected by a tag
  """
  
  def __init__(self,pdbfile) :
    self._vectors = []
    self._i = 0
    self._parse(open(pdbfile,'r').readlines())
    
  def __iter__(self) :
    return self
    
  def __len__(self) :
    return len(self._vectors)
    
  def next(self) :
    if self._i < len(self._vectors) :
      self._i += 1
      return self._vectors[self._i-1]
    else :
      raise StopIteration()
      
  def _parse(self,lines) :
    """
    Parse a list of PDB records
    
    Need to be overloaded by subclasses!
    """
    pass
    
  def _atomlines(self,lines) :
    """
    Generator that returns protein ATOM records of the backbone N atom,
    used by subclasses to iterate over only protein residues
    """
    i = 0
    while i < len(lines) :
      if lines[i][:6] == "ATOM  " and lines[i][12:16] == ' N  ' :
        yield lines[i]
      i += 1
      
  # Class variables and methods
  _vectorkeys = []
  _vectors = {}
   
  def vectorkeys(cls) :
    cls._init_vectors()
    return cls._vectorkeys
    
  def vectors(cls,tag) :
    cls._init_vectors()
    return cls._vectors[tag]
    
  def _init_vectors(cls) :
    """
    Initialises a dictionary of BondVectors sub classes and returns
    such a class based on a tag
    """
    if len(cls._vectorkeys) == 0 :
      def pred(c) :
        return inspect.isclass(c) and c.__module__ == cls.__module__ and issubclass(c,BondVectors) and c.__name__ != "BondVectors"
      for c in  inspect.getmembers(sys.modules[__name__],pred) :
        cls._vectorkeys.append(c[1].tag)
        cls._vectors[c[1].tag] = c[1]
        
  vectorkeys = classmethod(vectorkeys)      
  vectors = classmethod(vectors)
  _init_vectors = classmethod(_init_vectors)

  
class NHBonds(BondVectors) :
  """
  Parse bond vectors for backbone N-H groups
  Exclude proline residues
  """
  
  tag = "nh"

  def _parse(self,lines) :
    for s in self._atomlines(lines) :
      if s.find('PRO') > -1 : continue
      resnam = s[17:20]
      resid = s[22:26].strip()
      if int(resid) > 1 : self._vectors.append(BondVector(resnam,resid,'N','H'))
        
class MeBonds(BondVectors) :
  """
  Parse bond vectors for side-chain methyl groups
  """
  
  tag = "me"

  def _parse(self,lines) :

    for s in self._atomlines(lines) :
      resnam = s[17:20]
      resid = s[22:26].strip()
      if resnam == 'VAL' :
        self._vectors.append(BondVector(resnam,resid,"CB","CG1"))
        self._vectors.append(BondVector(resnam,resid,"CB","CG2"))
      elif resnam == 'THR' :
        self._vectors.append(BondVector(resnam,resid,"CB","CG2"))
      elif resnam == 'LEU' :
        self._vectors.append(BondVector(resnam,resid,"CG","CD1"))
        self._vectors.append(BondVector(resnam,resid,"CG","CD2"))
      elif resnam == 'ILE' :
        self._vectors.append(BondVector(resnam,resid,"CB","CG2"))
        self._vectors.append(BondVector(resnam,resid,"CG1","CD1"))
      elif resnam == 'MET' :
        self._vectors.append(BondVector(resnam,resid,"SD","CE"))
      elif resnam == 'ALA' :
        self._vectors.append(BondVector(resnam,resid,"CA","CB"))

class AlaBonds(BondVectors) :
  """
  Parse bond vectors for alanine side-chain methyl groups
  """
  
  tag = "ala"

  def _parse(self,lines) :
  
    for s in self._atomlines(lines) :
      resnam = s[17:20]
      resid = s[22:26].strip()
      if resnam == 'ALA' :
        self._vectors.append(BondVector(resnam,resid,"CA","CB"))

class AromaticBonds(BondVectors) :
  """
  Parse bond vectors for aromatic ring groups
  """
  
  tag = "ar"

  def _parse(self,lines) :
  
    for s in self._atomlines(lines) :
      resnam = s[17:20]
      resid = s[22:26].strip()
      if resnam == 'PHE' or resnam == 'TYR' :  
        self._vectors.append(BondVector(resnam,resid,"CD1","HD1")) 
        self._vectors.append(BondVector(resnam,resid,"CD2","HD2"))
        self._vectors.append(BondVector(resnam,resid,"CE1","HE1")) 
        self._vectors.append(BondVector(resnam,resid,"CE2","HE2"))  
      elif resnam == 'TRP' :  
        self._vectors.append(BondVector(resnam,resid,"CD1","HD1"))  
        self._vectors.append(BondVector(resnam,resid,"NE1","HE1")) 
        self._vectors.append(BondVector(resnam,resid,"CZ2","HZ2")) 
        self._vectors.append(BondVector(resnam,resid,"CE3","HE3")) 
        self._vectors.append(BondVector(resnam,resid,"CZ3","HZ3")) 
        self._vectors.append(BondVector(resnam,resid,"CH2","HH2")) 
      elif resnam == 'HID' :
        self._vectors.append(BondVector(resnam,resid,"CD2","HD2"))  
        self._vectors.append(BondVector(resnam,resid,"ND1","HD1"))  
        self._vectors.append(BondVector(resnam,resid,"CE1","HE1"))  
      elif resnam == 'HIE' :
        self._vectors.append(BondVector(resnam,resid,"CD2","HD2"))    
        self._vectors.append(BondVector(resnam,resid,"CE1","HE1")) 
        self._vectors.append(BondVector(resnam,resid,"NE2","HE2"))  
      elif resnam == 'HIP' :
        self._vectors.append(BondVector(resnam,resid,"CD2","HD2"))  
        self._vectors.append(BondVector(resnam,resid,"ND1","HD1"))  
        self._vectors.append(BondVector(resnam,resid,"CE1","HE1")) 
        self._vectors.append(BondVector(resnam,resid,"NE2","HE2"))  
          
          

class DictionaryBonds(BondVectors) :
  """
  Parse bond vectors for Bruschweilers dictionary
  """
  
  tag = "dic"

  def _parse(self,lines) :
  
    for s in self._atomlines(lines) :
      resnam = s[17:20]
      resid = s[22:26].strip()
      if resnam == 'VAL' :
        self._vectors.append(BondVector(resnam,resid,"CB","CG1"))
        self._vectors.append(BondVector(resnam,resid,"CB","CG2"))
      elif resnam == 'SER' :
        self._vectors.append(BondVector(resnam,resid,"CB","HB2"))
        self._vectors.append(BondVector(resnam,resid,"CB","HB3"))
      elif resnam == 'THR' :
        self._vectors.append(BondVector(resnam,resid,"CB","CG2"))
      elif resnam == 'ILE' :
        self._vectors.append(BondVector(resnam,resid,"CG1","CD1"))
      elif resnam == 'LEU' :
        self._vectors.append(BondVector(resnam,resid,"CG","CD1"))
        self._vectors.append(BondVector(resnam,resid,"CG","CD2"))
      elif resnam == 'MET' :
        self._vectors.append(BondVector(resnam,resid,"SD","CE"))
      elif resnam == 'ASN' :
        self._vectors.append(BondVector(resnam,resid,"ND2","HD21"))
        self._vectors.append(BondVector(resnam,resid,"ND2","HD22"))
      elif resnam == 'GLN' :
        self._vectors.append(BondVector(resnam,resid,"NE2","HE21"))
        self._vectors.append(BondVector(resnam,resid,"NE2","HE22"))
      elif resnam == 'PHE' :
        self._vectors.append(BondVector(resnam,resid,"CD1","HD1"))
      elif resnam == 'HID' or resnam == 'HIE' or resnam == 'HIP' :
        self._vectors.append(BondVector(resnam,resid,"CD2","HD2"))
      elif resnam == 'TYR' :
        self._vectors.append(BondVector(resnam,resid,"CD1","HD1"))
      elif resnam == 'PRO' :
        self._vectors.append(BondVector(resnam,resid,"CG","HG2"))
        self._vectors.append(BondVector(resnam,resid,"CG","HG3"))
      elif resnam == 'LYS' :
        self._vectors.append(BondVector(resnam,resid,"CB","HB2"))
        self._vectors.append(BondVector(resnam,resid,"CB","HB3"))
        self._vectors.append(BondVector(resnam,resid,"CG","HG2"))
        self._vectors.append(BondVector(resnam,resid,"CG","HG3"))
        self._vectors.append(BondVector(resnam,resid,"CD","HD2"))
        self._vectors.append(BondVector(resnam,resid,"CD","HD3"))
        self._vectors.append(BondVector(resnam,resid,"CE","HE2"))
        self._vectors.append(BondVector(resnam,resid,"CE","HE3"))
      elif resnam == 'ARG' :
        self._vectors.append(BondVector(resnam,resid,"CB","HB2"))
        self._vectors.append(BondVector(resnam,resid,"CB","HB3"))
        self._vectors.append(BondVector(resnam,resid,"CG","HG2"))
        self._vectors.append(BondVector(resnam,resid,"CG","HG3"))
        self._vectors.append(BondVector(resnam,resid,"CD","HD2"))
        self._vectors.append(BondVector(resnam,resid,"CD","HD3"))
        self._vectors.append(BondVector(resnam,resid,"NE","HE"))
      elif resnam == 'ASP' :
        self._vectors.append(BondVector(resnam,resid,"CB","HB2"))
        self._vectors.append(BondVector(resnam,resid,"CB","HB3"))
      elif resnam == 'GLU' :
        self._vectors.append(BondVector(resnam,resid,"CG","HG2"))
        self._vectors.append(BondVector(resnam,resid,"CG","HG3"))

