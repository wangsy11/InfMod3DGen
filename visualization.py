#Converting output modeling results into pdb format 


#import IMP
import string



for j in range(9):  # process in batch, file name can be modified
    k=j+1
    num='%d' %k
    str1='chr0'
    if (k==10):
        str1='chr'
    str2='visual'
    str3='.txt'
    str4='.pdb'
    
    file_object = open(str1+num+str3)
    try:
        t = file_object.readlines( )

    finally:
         file_object.close( )


    file_object = open(str2+num+str4, 'w')

    for i in range(len(t)):
     Enter = '\n'
     Pos = t[i].index(Enter)
     file_object.writelines('ATOM     ')
     if (i<9):
       file_object.writelines(' ')
     file_object.writelines(str(i+1))
     file_object.writelines(' C    LIG A           ')
     file_object.writelines(t[i][0:6])
     file_object.writelines('  ')
     file_object.writelines(t[i][6:12])
     file_object.writelines('  ')
     file_object.writelines(t[i][12:17])
     file_object.writelines('  1.00 75.00    \n') 

    for i in range(len(t)-1):
     file_object.writelines('CONECT   ')
     if (i<9):
       file_object.writelines(' ')
     file_object.writelines(str(i+1))  
     file_object.writelines('  ')
     if (i+1<9):
       file_object.writelines(' ')
     file_object.writelines(str(i+2)) 
     file_object.writelines('\n')

    file_object.writelines('END\n') 

    file_object.close( )
           


