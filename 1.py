#Script for creating solvent boxes

#Import method for interacting with the operating system via command line
import os 

#Initializing Lists
atomTypes = []
bondTypes = []
angleTypes = []
dihedTypes = []

molecData = []
molecBonds = []
molecAngles = []
molecDiheds = []

atomList = []
bondList = []
angleList = []
dihedList = []

#Beginning of Input Section.........................................................................................
outputFileName = 'Solvent.data' 	#Self-explanatory
solvDens = 			#Solvent Density: Solvent Molecule Number / Angstrom^3     
BoxDimension = 25.0                   	#Length Dimension of an equilibrated cubic box in Angstroms
ExpansionFactor =     			#Multiplicative factor for BoxDimension that initializes the molecules in an expanded box

#One example of each type in the following sections:

#Atom types ([index, atomic mass])
atomTypes.append([1, 15.999])

#Bond types ([index, k, r0])
bondTypes.append([1, 600, 1.5])

#Angle types ([index, k, theta0])
angleTypes.append([1, 100, 109])

#Dihedral types ([index, V1, V2, V3])
dihedTypes.append([1, 0.0, 1.0, 1.0])

#Molecule Atom List ([Index, 1, Atom Type, q, x, y, z])
molecData.append([1, 1, 1, 0.5, -2.76, 3.79, 9.82])

#Molecule Bond List ([Index, Bond Type, Atom Index 1, Atom Index 2])
molecBonds.append([1, 1, 1, 2])

#Molecule Angle List ([Index, Angle Type, Atom Index 1, Atom Index 2, Atom Index 3])
molecAngles.append([1, 1, 1, 2, 3])

#Molecule Dihedral List ([Index, Angle Type, Atom Index 1, Atom Index 2, Atom Index 3])
molecDiheds.append([1, 1, 1, 2, 3, 6])

#End of inputs........................................................................................

#Determining length of lists
atomNum = len(molecData)
bondNum = len(molecBonds)
angleNum = len(molecAngles)
dihedNum = len(molecDiheds)

totSolv = int((BoxDimension**3) * (solvDens) ) #Determines number of solvent molecules
boxSide = BoxDimension * ExpansionFactor
dims = [-boxSide/2.0, boxSide/2.0, -boxSide/2.0, boxSide/2.0, -boxSide/2.0, boxSide/2.0 ] #Box dimensions = [x0, x1, y0, y1, z0, z1]

#Creating total list of solvent molecules
for i in range(totSolv):
	for j in range(atomNum):
		atomList.append( [ molecData[j][0]+i*atomNum, molecData[j][1]+i, molecData[j][2], molecData[j][3], molecData[j][4], molecData[j][5], molecData[j][6] ] )
	for j in range(bondNum):
		bondList.append( [ molecBonds[j][0] + i*bondNum, molecBonds[j][1], molecBonds[j][2]+i*atomNum, molecBonds[j][3]+i*atomNum ] )
	for j in range(angleNum):
		angleList.append( [ molecAngles[j][0] + i*angleNum, molecAngles[j][1], molecAngles[j][2]+i*atomNum, molecAngles[j][3]+i*atomNum, molecAngles[j][4]+i*atomNum ] )
	for j in range(dihedNum):
		dihedList.append( [ molecDiheds[j][0] + i*dihedNum, molecDiheds[j][1], molecDiheds[j][2]+i*atomNum, molecDiheds[j][3]+i*atomNum, molecDiheds[j][4]+i*atomNum, + molecDiheds[j][5]+i*atomNum ] )

#Deleting files that pertain to Packmol
os.system('rm tempCoords.xyz')
os.system('rm temp.inp')
os.system('rm packCoords.xyz')

#Creating temporary molecule coordinates for Packmol (tempCoords.xyz)
solvFile = open('tempCoords.xyz', 'a')
solvFile.write(str(atomNum) + '\n')
solvFile.write('Atoms\n')
for i in range(atomNum):
	solvFile.write('%s	%s	%s	%s\n' % (molecData[i][2], molecData[i][4], molecData[i][5], molecData[i][6]) )
solvFile.close()

#Creating Packmol script (temp.inp)
packmol = open('temp.inp', 'a')
packmol.write('tolerance 2.0\n')
packmol.write('output packCoords.xyz\n')
packmol.write('filetype xyz\n')
packmol.write('structure tempCoords.xyz\n')
packmol.write('	number ' + str(totSolv) + '\n')
packmol.write('	inside box %s %s %s %s %s %s\n' % (dims[0], dims[2], dims[4], dims[1], dims[3], dims[5]))
packmol.write('end structure\n')
packmol.close()

os.system('/fs/home/ns728/bin/packmol < temp.inp')

#Collecting Packmol Coordinates
packedCoords = open('packCoords.xyz')
inStr = packedCoords.readline()
inStr = packedCoords.readline()
for i in range(atomNum*totSolv):
	inStr = packedCoords.readline()
	inList = inStr.split()
	atomList[i][4] = float(inList[1])
	atomList[i][5] = float(inList[2])
	atomList[i][6] = float(inList[3])

#Creating Data File
os.system('rm ' + outputFileName)

outFile = open(outputFileName, 'a')

outFile.write('LAMMPS Description \n\n')

outFile.write('     %s atoms\n' % (atomNum*totSolv))                         
outFile.write('     %s bonds\n' % (bondNum*totSolv))
outFile.write('     %s angles\n' % (angleNum*totSolv))
outFile.write('     %s dihedrals\n\n' % (dihedNum*totSolv))

outFile.write('     %s atom types\n' % (len(atomTypes)))
outFile.write('     %s bond types\n' % (len(bondTypes)))
outFile.write('     %s angle types\n' % (len(angleTypes)))
outFile.write('     %s dihedral types\n\n' % (len(dihedTypes)))


#Print dimensions
outFile.write('%s %s xlo xhi\n' % (dims[0], dims[1]))
outFile.write('%s %s ylo yhi\n' % (dims[2], dims[3]))
outFile.write('%s %s zlo zhi\n\n' % (dims[4], dims[5])) 

#Printing Types
outFile.write('Masses\n\n')  
for i in range(len(atomTypes)):
	outFile.write(' %s %s\n' % (atomTypes[i][0], atomTypes[i][1]) )
outFile.write('\n')

outFile.write('Bond Coeffs\n\n')  
for i in range(len(bondTypes)):
	outFile.write(' %s %s %s\n' % (bondTypes[i][0], bondTypes[i][1], bondTypes[i][2]) )
outFile.write('\n')

outFile.write('Angle Coeffs\n\n')  
for i in range(len(angleTypes)):
	outFile.write(' %s %s %s 0 0\n' % (angleTypes[i][0], angleTypes[i][1], angleTypes[i][2]) )
outFile.write('\n')

outFile.write('Dihedral Coeffs\n\n')  
for i in range(len(dihedTypes)):
	outFile.write(' %s %s %s %s 0.0\n' % (dihedTypes[i][0], dihedTypes[i][1], dihedTypes[i][2], dihedTypes[i][3]) )
outFile.write('\n')

#Printing Atomic Coordinates
outFile.write('Atoms\n\n')
for i in range(atomNum*totSolv): #edited for charge
	outFile.write('  %s	%s	%s	%s	%s	%s	%s\n' % (atomList[i][0], atomList[i][1], atomList[i][2], atomList[i][3], atomList[i][4], atomList[i][5], atomList[i][6]))   
outFile.write('\n')

#Printing Bond, Angle, and Dihedral Data
outFile.write('Bonds\n\n')
for i in range(bondNum*totSolv):
	outFile.write('  %s	%s	%s	%s\n' % (bondList[i][0], bondList[i][1], bondList[i][2], bondList[i][3]))
outFile.write('\n')

outFile.write('Angles\n\n')
for i in range(angleNum*totSolv):
	outFile.write('  %s	%s	%s	%s	%s\n' % (angleList[i][0], angleList[i][1], angleList[i][2], angleList[i][3], angleList[i][4]))
outFile.write('\n')

outFile.write('Dihedrals\n\n')
for i in range(dihedNum*totSolv):
	outFile.write('  %s	%s	%s	%s	%s	%s\n' % (dihedList[i][0], dihedList[i][1], dihedList[i][2], dihedList[i][3], dihedList[i][4], dihedList[i][5]))

outFile.close()
#Cleaning Packmol Information
os.system('rm tempCoords.xyz')
os.system('rm temp.inp')
os.system('rm packCoords.xyz')



