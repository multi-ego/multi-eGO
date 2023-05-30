from ast import arg
from inspect import Parameter
import os
import subprocess
import pandas as pd
import sys, getopt
import shlex
import shutil 
import ntpath
import time
import random


REFERENCE='./inputs'
#MD='./MD'

aforismi=[
    '“Imparerai a tue spese che nel lungo tragitto della vita incontrerai tante maschere e pochi volti.”',
    "“Mi si fissò invece il pensiero ch'io non ero per gli altri quel che finora, dentro di me, m'ero figurato d'essere.”",
    "“È molto più facileessere un eroe che un galantuomo. Eroi si può essere una volta tanto; galantuomini, si dev'esser sempre.”",
    "“Basta che lei si metta a gridare in faccia a tutti la verità. Nessuno ci crede, e tutti la prendono per pazza!”",
    "“La civiltà vuole che si auguri il buon giorno a uno che volentieri si manderebbe al diavolo; ed essere bene educati vuol dire appunto esser commedianti.”",
    "“Sapete che cosa significa amare l'umanità? Significa soltanto questo: essere contenti di noi stessi. Quando uno è contento di sé stesso, ama l'umanità.”",
    "“La facoltà d'illuderci che la realtà d'oggi sia la sola vera, se da un canto ci sostiene, dall'altro ci precipita in un vuoto senza fine, perché la realtà d'oggi é destinata a scoprire l'illusione domani. E la vita non conclude. Non può concludere. Se domani conclude, è finita.”",
    "“Ciò che conosciamo di noi è solamente una parte, e forse piccolissima, di ciò che siamo a nostra insaputa.”",
    "“Confidarsi con qualcuno, questo sì è veramente da pazzo!”",
    "“Se noi conosciamo che errare è dell'uomo non è crudeltà sovrumana la giustizia?”",
    "“Io non l'ho più questo bisogno, perché muojo ogni attimo io, e rinasco nuovo e senza ricordi: vivo e intero, non più in me, ma in ogni cosa fuori.”",
    "“Notiamo facilmente i difetti altrui e non ci accorgiamo dei nostri.”",
    "“Ogni cosa finché dura porta con sé la pena della sua forma, la pena d'esser così e di non poter essere più altrimenti.”",
]

def print_message(message):
    print('----------------------------------------------------')
    print('')
    print('LUIGI Message:)   '+message)
    print('')
    print('----------------------------------------------------')

def start_message(message):
    print('')
    print('LUIGI Message:)   '+message)
    print('----------------------------------------------------')

def print_inside_message(message):
    print('')
    print(message)


def end_message():
    print('')
    print('----------------------------------------------------')



def print_message_help():
    print('')
    print('LUIGI Helps the needy:)')
    print('----------------------------------------------------')
    print('')
    print('')
    print('luigi.py usage:')
    print('     --protein or -p (optional) [protein_name or PATH/protein_name]: give the protein_name or path of your protein_name.pdb file')
    print('     --gmx or -g (optional) [/PATH_TO_GMX/gmx or /PATH_TO_GMX/gmx_mpi]: If gmx path (particular version) was not specified, luigi will search and use the first found')
    print('     --remove or -r (optional) [HOH, BNZ,...]: molecules to be removed from the pdb file')
    print('     --md (optional)[md_name or PATH/md_name]: path where to find md simulation')
    print('     --md_name (optional): The new name of the md file (and the directory where it will be contained). If not specified luigi will use the name given in --md')
    print('     --input_name (optional): This will tell Luigi where to copy the md ensemble files in case --protein was not specified (because the reference structure was already made)')
    print('     --rc (optional)  random coil of reference structure')
    print('') 
    print('')
    print('----------------------------------------------------')
    print('')


def print_warning(message):
    print('----------------------------------------------------')
    print('')
    print('LUIGI is worried :S    '+message)
    print('') 
    print('----------------------------------------------------')

def print_error(error):
    print('----------------------------------------------------')
    print('')   
    print('Luigi found an ERROR: '+error)
    print('')
    print('----------------------------------------------------')

def print_goodbye():
    print('-----------------------------===========+===+++++++++++++********#########################')
    print('-----------------------------=========++++**##########**********##########################')
    print('-----------------------------=-==+**#%%%%####*****#####%%%@%#########%%%%%################')
    print('----------------------=---====*#%@@@%%%%##**#******+****##%%@@@@@%%%%%%%##################')
    print('-------------------------==+#%@@@%%%###***+*+++=+*+*+++***##%@@@@@@@@@@%#%################')
    print('-----------------------==+%@@%%%###***+++========++=+++++**####%%@@@@@@@@%################')
    print('---------------------=+#%@@%###****+++++=============++=++++***###%@@@@@@@%###############')
    print('--------------------=*%%@@%###***+++====================+++++****###%@@@@@@%##############')
    print('------------------=*%@@@@%#***+++==========--============++++******##%%@@@@@@%############')
    print('-----------------=#@@@@@%#***++========-===--===--========+++++++****##%@@@@@@%#####%%%%%#')
    print('----------------+@@@@@%##**+++=====-==-=------------=======++++++++***##%%@@@@@###%%%%%%%%')
    print('---------------+@@@@@%#*+++=+=======--------==----------=========+**+*#%%%@@@@@%%%%%%%%%%%')
    print('--------------=%@@@@@##**++===-=---=--------------------=====++==+++***#%@@@@@@@%%%%@@@@@@')
    print('--------------+@@@@@%##*+++===----------------------------====+++++++*#%%%%@@@@@%%%%@@@@@@')
    print('-------------=#@@@@%%##+++===--------------------------------==+++++***#%%##@@@@%%%%%%@@@@')
    print('-------------=%@@@@%%#*+++=--------------------------------=====++++***#%%%%#@@@%%%%%%%%%%')
    print('-------------+@@@@@@###*++==-------------------------------====++++****#%%#%%%@@%#%%%%%%%%')
    print('------------=%@@@@@%###++===-=--------------------------======+++++****#%%%##%@@%######%%%')
    print('------------+@@@@@%%%#*++=====-==------------------=========++=*+++***####@@%%@@#########%')
    print('------------+@@@@@%%%**+======-=--------------------=======++++**********#%%%%@@########%%')
    print('------------+@@@@%%%##+++========------------------=-======+++++***+**###**###@@###*######')
    print('------------=%@@@%###*+=====------------------------======+++++++**++#*====+*#+=====+#####')
    print('------------=*@@@@#**+++==-===-------------=--------======++===+#+===+==----=+*+=---=*####')
    print('-------------=@@@@@@%*+====--=-------------------========+====+*=---==+*#*=--=+++=--=+%###')
    print('-------------+@@@@@@@@%#*=====---===----------------===++=====+==---#+=+*@#=--=+*+=--=#%##')
    print('-------------#@@@@@@@@@@@%#+======+======-========+*##%#+++=--===-=**#+=+======+++*--=%@##')
    print('------------=%@@@@@@@@@@@@@@%++===+=====++*#%%%@@@@%#*+=+%+=-------=#%+=====+=-===++-+#@@@')
    print('------------=@@@@@#@##@@@@@@@@#+=+*+==+*#%@@@%##*+=+**++#%+==-------=-==---==++===-=++##@@')
    print('------------=@%##%%**+###*#@@@@+=====+#%%%#%%@%%%%@@@@#*%++%#++---------=======*+==-+##%@@')
    print('------------+@##**#%%##**+++#@@+=--==++#%%%%%**@@@@%%%@@#**%@%+==---------===--=*+==+###@@')
    print('------------*@%###********++#@%==--===++=+*+==+++++=+++#%==+@#===---------===--=+++=+%%#@@')
    print('-----------=#@@##%###**###*#%%*=----=====*%%%#*++=====++@+==#*+---------==++===**#*+%@%%@@')
    print('-------===+#@@@%#*###**++##@@%+---====+--=+#%%%%##*#*++=#%+==---------==+***++*#%*##%@@@@@')
    print('---===*#%@@@@@@%*++======*%@@*=---=====-----==+*#*+=====+*=---------==+**+*+=+%%+**@@@@@@@')
    print('==+#%@@@@@@@@@@%*+======+*%@@*=---=====-----=======-===+===---------=++++++=+%#**#@@@@@@@@')
    print('#@@@@@@@@@@@@@@@#++=++*#%@%@@+=---====-==-=-====----=+*+==---------=====++=+#+*+#@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@#**##%%@#*@*==---=========+++====--++===--------=====++++=#+=+*@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@#*#@###==----====--==+=**+====+*===---------====+====*==*%@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@%+#@@@%+==---====--=+====*#*+++%+==---------==-===-=+#=+#@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@**+=+@@@@+=========+##+====+*%##@#+==------------==-=++=+#@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@%%#**%@@%*+*++#%@%#+=-=====++#@@#+==------------=-=*==-+#@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@##@@@@@@@@@%*==-------=-===%@#+===------------=+==-=+@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*===-=---======-+@@@#+=-------------=====+*%@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%%##+***%%*%%%%@@@@*=--------------==+=++@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%##**##*#+**+**+++*@@@@@@+--------------==+**%@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@%@@%#*#****#*+++*+=====+%@@@@@#=---------------=+**@%@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@*++#@%===-===============*@@@@@@+----------------==---*@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@*+=++===------========--=#%@@@@#=---------------------+@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@%=--=+==-------=---------*@@@@@%=----------------------*@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@*=-=------------------=+%@@@@@@=---------------------=*@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@#+=------------=====+%@@@@@@@@=------------------=*#@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@*=-----------==+*%@@@@@%#*++=----------------+#@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%=--------=*%%@@@@%#+==-----------------=+*%@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%+==--=+*%@@@@@@@=-------------------=*%@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%##%@@@@@@@@@@*----------------=+#@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@+==-----==---=*#@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###*=---**+*%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%%%@@@@@@@@@@@@@@@@@#+*%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('')
    print('')
    print('')
    print('----------------------------------------------------')
    print('')
    print(':-)  LUIGI says goodbye  (-:')
    print('')
    print(aforismi[random.choice(range(len(aforismi)))])
    print()
    print('')
    print('----------------------------------------------------')
    print('')
    print('')
    print('')

def print_hello():
    print('')
    print('')
    print('')
    print('----------------------------------------------------')
    print('')
    print(':-)  LUIGI Multi-identity menager  (-: ')
    print('')
    print('----------------------------------------------------')
    print('')
    print('')
    print('')
    time.sleep(1)



def check_pdb(protein):

    """Check whther protein file is a pdb file: finishes with '.pdb' """

    size=len(protein)
    if protein[size-4:size]=='.pdb':
        return True
    else:
        return False

def is_tool(name):

    """Check whether `name` is on PATH and marked as executable."""

    return shutil.which(name) is not None


def read_options(opts, parameters):

    """Read the options given to luigi and puts them in parameters"""

    for opt, arg in opts:
        if opt in ("--help", "-h"):
            print_message_help()
            sys.exit()

        elif opt in ("--rc"):
            print_message('Luigi will calculate Random Coil')
            parameters['rc']='rc'


        elif opt in ("--protein", "-p"):
            if not arg:
                print_message('Provide a protein name (or path if in a different directory) of your protein_name.pdb file: --protein protein_name')
                sys.exit()

            elif arg:
                prot=arg
                tail, head=ntpath.split(arg)
                head_split=head.split('.')
                dim=len(head_split)
                #controls if the path+name inserted is correct and accessible
                if os.path.exists(prot)==False:
                    print_error("Protein file 2: "+prot+" cannot be found or doesn't exists")
                    sys.exit()
                else:
                    #controls if the name finishes with .pdb, otherwise Luigi adds it
                    if dim==1:
                        prot=prot+'.pdb'
                        print_warning('missing ".pdb", Luigi will add it')
                    elif dim==2 and head_split[1]=='pdb':
                        print_message('correctly inserted pdb file')
                    else:
                        print_error('Something went wrong when selecting protein file. Please insert "protein name" or "protein_name.pdb"')
                        sys.exit()


                    tail, head=ntpath.split(prot)
                    if tail=='':
                        parameters['protein_file']=head
                        parameters['protein_name']=head.split('.')[0]
                        parameters['protein_path']='.'
                    else:
                        parameters['protein_file']=head
                        parameters['protein_name']=head.split('.')[0]                        
                        parameters['protein_path']=tail   

            else:
                print_error('Something horrible happend in finding your protein')
                sys.exit()
            print_message('PROTEIN  IS: \n  \n protein file='+parameters['protein_file']+'\n protein name='+ parameters['protein_name']+'\n protein path='+parameters['protein_path'])
                    
        elif opt in ("--gmx", "-g"):
            start_message('Luigi is going to define the gmx command')
            if  arg:
                if os.path.exists(arg)==True:
                    parameters['gmx_path']=arg
                else:
                    print_warning("The given gmx command path was not found. Luigi will search for gmx command and use that ")
                    if is_tool('gmx'):
                        print_inside_message('luigi is going to search for gmx command')
                        parameters['gmx_path']=shutil.which('gmx')
                        print_inside_message('Luigi found a command. Is going to use '+ parameters['gmx_path'])

                    elif is_tool('gmx_mpi'):
                        print_inside_message('luigi is going to search for gmx_mpi command')
                        parameters['gmx_path']=shutil.which('gmx_mpi')
                        print_inside_message('Luigi found a command. Is going to use '+ parameters['gmx_path'])
                  
            else:
                if is_tool('gmx'):
                    print_warning('No gmx command was selected. Luigi is going to search for gmx command')
                    parameters['gmx_path']=shutil.which('gmx')
                    print_inside_message('Luigi found a command. Is going to use '+ parameters['gmx_path'])

                elif is_tool('gmx_mpi'):
                    print_warning('luigi is going to search for gmx_mpi command')
                    parameters['gmx_path']=shutil.which('gmx_mpi')
                    print_inside_message('Luigi found a command. Is going to use '+ parameters['gmx_path'])
                else:
                    print_error("gmx command not found.Please install gmx ")
                    sys.exit()
            
            end_message()

        elif opt in ("--remove", "-r"):

            """read the molecules to be removed in the format HOH,BNZ and convert it in a usable string 'HOH\|BNZ\|...'"""

            if arg:
                list_arg="["+arg+"]"
                list_appo=list_arg[1:-1].split(',')
                parameters['remove']=list_appo
                to_remove_list="'"+parameters['remove'][0]
                for i in range(1,len(parameters['remove'])):
                    to_remove_list=to_remove_list+"\|"+parameters['remove'][i]
                to_remove_list+="'"
                parameters['to_remove_list']=to_remove_list

                print_message(f'The molecules that will be removed from the pdb file are: {list_appo}')
            else:
                print_error('Provide a series of names to be removed "HOH,BNZ,..."')
                sys.exit()

        elif opt in ("--md"):
            if not arg:
                print_message('Provide a md.gro file (with path) of your md_name.gro file: --md_path md.gro')
                sys.exit()

            elif arg:
                md=arg
                tail, head=ntpath.split(arg)
                head_split=head.split('.')
                dim=len(head_split)

                #controls if the path+name inserted is correct and accessible
                if os.path.exists(md)==False:
                    print_error("Protein file 2: "+md+" cannot be found or doesn't exists")
                    sys.exit()

                else:
                    #controls if the name finishes with .pdb, otherwise Luigi adds it
                    if dim==1:
                        md=md+'.gro'
                        print_warning('missing ".gro", Luigi will add it')
                    elif dim==2 and head_split[1]=='gro':
                        print_message('correctly inserted gro file for md simulation')
                    else:
                        print_error('Something went wrong when selecting protein file. Please insert "md_name" or "md_name.gro"')
                        sys.exit()


                    tail, head=ntpath.split(md)
                    if tail=='':
                        parameters['md_name']=head.split('.')[0]
                        parameters['md_path']='.'
                    else:
                        parameters['md_name']=head.split('.')[0]                        
                        parameters['md_path']=tail           

                    parameters['md_new_name']=parameters['md_name']
            else:
                print_error('Something horrible happend in finding your md simulation')
                sys.exit()

            print_message('MD  IS: \n  \n md file='+parameters['md_name']+'.gro \n md name='+ parameters['md_name']+'\n md path='+parameters['md_path'])
            
            
            
        elif opt in ("--md_name"):
            if not arg:
                print_message('Provide a name for the md: --md_name name')
                sys.exit()

            elif arg:
                parameters['md_new_name']=arg
            
            print_message('Md name will not be the default one but instead: \n \n'+parameters['md_new_name'])


        elif opt in ("--input_name"):
            if not arg:
                print_message('Provide a name for the input directory. If specified it must be equal to protein_name: --p_name name')
                sys.exit()

            elif arg:
                parameters['input_name']=arg
        

            ###END of for###

 

    #prepare parameters of imput directory names

    #if protein name and input name are specified control that they are the same
    if parameters['protein_name']!='none' and parameters['input_name']!='none':

        if parameters['protein_name']!=parameters['input_name']:
            print_error('If protein is specified, the input directory is not necessary. However, if both are selected they must be the same ')
            sys.exit()
        else:
            parameters['ref_path']          =REFERENCE+'/'+parameters['protein_name']+'/reference' 
            parameters['protein_input_path']=REFERENCE+'/'+parameters['protein_name']
            parameters['md_input_path']     =REFERENCE+'/'+parameters['protein_name']+'/md_ensembles' 

    #if only protein name is specified
    elif parameters['protein_name']!='none' and parameters['input_name']=='none':

        parameters['ref_path']          =REFERENCE+'/'+parameters['protein_name']+'/reference' 
        parameters['protein_input_path']=REFERENCE+'/'+parameters['protein_name']
        parameters['md_input_path']     =REFERENCE+'/'+parameters['protein_name']+'/md_ensembles' 

     #if protein name is not specified and input is
    elif parameters['protein_name']=='none' and parameters['input_name']!='none':
        if os.path.exists(REFERENCE+'/'+parameters['input_name']):
            parameters['ref_path']          =REFERENCE+'/'+parameters['input_name']+'/reference' 
            parameters['protein_input_path']=REFERENCE+'/'+parameters['input_name']
            parameters['md_input_path']     =REFERENCE+'/'+parameters['input_name']+'/md_ensembles' 
        else:
            print_error(f"{REFERENCE}/{parameters['input_name']} was not found")
            sys.exit()

    else:
        print_error('Option not taken into account by luigi')
        sys.exit()


    #search for the gmx_command if was not specified
    if parameters['gmx_path']=='none':
        if is_tool('gmx'):
            print_warning('No gmx was selected. Luigi is going to search for gmx command')
            parameters['gmx_path']=shutil.which('gmx')
            print_message('Luigi found a command. Is going to use '+ parameters['gmx_path'])

        elif is_tool('gmx_mpi'):
            print_warning('luigi is going to search for gmx_mpi command')
            parameters['gmx_path']=shutil.which('gmx_mpi')
            print_message('Luigi found a command. Is going to use '+ parameters['gmx_path'])
        else:
            print_error("gmx command not found.Please install gmx ")
            sys.exit()
            
    return parameters


def Make_REF_directory(parameters):

    """Generate REFERENCE directory if does not exists in this directory"""

    if os.path.exists(REFERENCE)==False:

        os.mkdir(REFERENCE)


    if os.path.exists(REFERENCE+'/'+parameters['protein_name'])==False:

        os.mkdir(REFERENCE+'/'+parameters['protein_name'])


    if os.path.exists(parameters['ref_path'])==False:

        os.mkdir(parameters['ref_path'])

    shutil.copytree("multi-ego-basic.ff", f"{parameters['ref_path']}/multi-ego-basic.ff", dirs_exist_ok=True)

def Make_MD_directory(parameters):

    """Generate MD directory if does not exists in this directory"""

    if os.path.exists(parameters['protein_input_path'])==False:

        print_error(f"{parameters['protein_input_path']} doesn't exists. First make specifying to luigi a reference structure with --protein  " )
        sys.exit()   


    elif parameters['protein_input_path']+'/md_ensembles'==parameters['md_input_path']:
        if os.path.exists(parameters['md_input_path'])==False:
            os.mkdir(parameters['md_input_path'])
            print_message(f"Luigi created {parameters['md_input_path']}")

        else:
            print_message(f"Directory {parameters['md_input_path']} already exists")

    else:
        print_error(f"{parameters['protein_input_path']}/md_ensembles != {parameters['md_input_path']}")
        sys.exit()

    

def Copy_protein(parameters):

    '''Copies the protein from the directory it is in to input/protein/reference'''

    if os.path.exists(parameters['ref_path'])==True:
        
        source=   f"{parameters['protein_path']}/{parameters['protein_file']}"
        destination=  f"{parameters['ref_path']}/{parameters['protein_file']}"
        shutil.copy(source, destination)
        print_message("File copied successfully.")

    else:
        print_error(f"directory {parameters['ref_path']} does not exists")   
        sys.exit()  


def Copy_md(parameters):

    '''Copies the md files from the directory it is in input/protein/md_ensemble/md_new_name'''

    if os.path.exists(parameters['md_input_path'])==True:

        destination_dir=f"{parameters['md_input_path']}/{parameters['md_new_name']}"

        #if the path already exists firstly it will be removed and then remade
        if os.path.exists(destination_dir): 
            shutil.rmtree(destination_dir)
            print_message(f'Luigi succesfully removed {destination_dir}')

        os.mkdir(destination_dir)
        print_message(f'Luigi succesfully made {destination_dir}')


        #Copy md data in inputs/protein/md_ensembles directory
        source_base=f"{parameters['md_path']}/{parameters['md_name']}"          #define a base name old_PATH/md_name     and later add '.gro', '.xtc', ecc...
        destination_base=destination_dir+'/'+parameters['md_new_name']          #define a base name new_PATH/md_new_name and later add '.gro', '.xtc', ecc...

        start_message(f"Copy all files from {parameters['md_path']} to {destination_dir}")

        #.gro
        source=source_base+'.gro'
        destination=destination_base+'.gro'
        shutil.copyfile(source,destination)
        print_inside_message(f'Luigi succesfully copied {source} in {destination}')
        
        #.tpr
        source=source_base+'.tpr'
        destination=destination_base+'.tpr'
        shutil.copyfile(source,destination)
        print_inside_message(f'Luigi succesfully copied {source} in {destination}')

        #.xtc
        source=source_base+'.xtc'
        destination=destination_base+'.xtc'
        shutil.copyfile(source,destination)
        print_inside_message(f'Luigi succesfully copied {source} in {destination}')

        #.top
        source=parameters['md_path']+"/topol.top"
        destination=destination_dir +"/topol.top"
        shutil.copyfile(source,destination)
        print_inside_message(f'Luigi succesfully copied {source} in {destination}')


        #copy force field in topology
        with open(parameters['md_path']+"/topol.top", 'r') as fp:
            lines = fp.readlines()
            target_line='none'
            for row in lines:
                word = '.ff'
                if row.find(word) != -1:
                    target_line=row
                    break

        if target_line=='none':
            print_error('No force field was found in '+parameters['md_path']+"/topol.top. Please control the topology file")
            sys.exit()
        else:
            ff_appo=target_line.split(' ')[1]
            ff, altro=ntpath.split(ff_appo)
            ff_path,ff_name=ntpath.split(ff) 
            ff_path=ff_path[1:]

            #if force field is a default one
            if ff_path=='':
                print_message('The force field is one of the standards of gromacs. Luigi will not copy it in the MD directory')

            #if the force field is in the working directory
            elif ff_path=='.':

                source=parameters['md_path']+'/'+ff_name
                destination=destination_dir+'/'+ff_name
                shutil.copytree(source, destination)
                print_inside_message(f'Luigi succesfully copied {source} in {destination}')

            #if the force field is in the previous directory
            elif ff_path=='..':

                source=parameters['md_path']+'/'+ff_path+'/'+ff_name
                destination=destination_dir+'/'+ff_name
                shutil.copytree(source, destination)
                print_inside_message(f'Luigi succesfully copied {source} in {destination}')

            #use the full path
            else:
                source=ff
                destination=destination_dir+'/'+ff_name
                shutil.copytree(source, destination)
                print_inside_message(f'Luigi succesfully copied {source} in {destination}')

    else:
        print_error(f"Directory {parameters['md_input_path']} doesn't exists")
        sys.exit()

def Remove_molecules(parameters):

    '''Removes unnecessary moleclues from pdb and rename protein.pdb to protein_clean.pdb'''

    if parameters['to_remove_list']!='none':
        try:
            os.system(f"grep -v {parameters['to_remove_list']} {parameters['ref_path']}/{parameters['protein_file']} > {parameters['ref_path']}/{parameters['protein_name']}_clean.pdb ")
        except:
            print_error("Cannot remove unnecessary molecules")
    else:
        try:
            os.system('cp '+parameters['ref_path']+'/'+parameters['protein_name']+'.pdb '+parameters['ref_path']+'/'+parameters['protein_name']+'_clean.pdb')
        except:
            print_error(f"something went wrong while making {parameters['protein_name']}_clean.pdb")
               

def Make_topology(parameters):

    '''Use pdb2gmx to create topology with defined force field and water model. Then copy posre.tip, topol.top and .gro files to REFERENCE'''

    start_message(f'Luigi is making the topology of reference structure: {parameters["protein_name"]}')

    #make pdb2gmx with default 1 and 1 for force field and water model
    try:
        entry="echo 1; echo 1"    #Choosing default parameters 1 and 1 for the forcefield
        p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True)
        subprocess.run(shlex.split(f"{parameters['gmx_path']} pdb2gmx -f {parameters['ref_path']}/{parameters['protein_name']}_clean.pdb -o {parameters['ref_path']}/{parameters['protein_name']}_clean.gro -ignh"), check = True, stdin=p1.stdout)
    except FileNotFoundError as exc:
        print_error(f"Process failed because the executable could not be found.\n{exc}")
    except subprocess.CalledProcessError as exc:
        print_error("pdb2gmx doesn't work")
        sys.exit()
    
    #move .top to reference
    source='./topol.top'
    destination=parameters['ref_path']+'/topol.top'
    shutil.move(source, destination)
    print_inside_message(f'Luigi succesfully moved {source} in {destination}')

    #move posre.itp to reference
    source='./posre.itp'
    destination=parameters['ref_path']+'/posre.itp'
    shutil.move(source, destination)
    print_inside_message(f'Luigi succesfully moved {source} in {destination}')
    
    #make the correct .gro file for multi-ego (naming and enumeration in the starting pdb is not correct)
    try:
        subprocess.run(shlex.split(f"{parameters['gmx_path']} editconf -f {parameters['ref_path']}/{parameters['protein_name']}_clean.gro -o {parameters['ref_path']}/{parameters['protein_name']}_ordered.pdb"), check = True)
        print_inside_message('Luigi made the correct .pdb file')
    except FileNotFoundError as exc:
        print_error(f"Process failed because the executable could not be found.\n{exc}")
    except subprocess.CalledProcessError as exc:
        print_error("editconf doesn't work")
        sys.exit()

    #change the pdb file in reference directory with the new and corrected one
    old_name=f"{parameters['ref_path']}/{parameters['protein_name']}_ordered.pdb"
    new_name=f"{parameters['ref_path']}/{parameters['protein_name']}.pdb"
    os.rename(old_name, new_name)
    print_inside_message('Luigi moved new.pdb file to old.pdb file in reference directory')

    #change name of cleaned .gro file
    old_name=f"{parameters['ref_path']}/{parameters['protein_name']}_clean.gro"
    new_name=f"{parameters['ref_path']}/{parameters['protein_name']}.gro"
    os.rename(old_name, new_name)
    print_inside_message(f'Luigi moved {old_name} file to .gro file in reference directory')

   
    #remove protein_cleaned.pdb
    file_to_remove=f"{parameters['ref_path']}/{parameters['protein_name']}_clean.pdb"
    os.remove(file_to_remove)
    print_inside_message(f'Luigi removed {file_to_remove}')

    end_message()


def Make_xtc_noH(parameters):

    """
    Prepare all the necessary to create an xtc file of inly the protein without hydrogens. In order:
        - convert the starting tpr to a tpr file only with protein
        - convert the starting tpr to a tpr file only with protein without hydrogens
        - from the tpr without hydrogen make a gro file (needed for MDmat)
        - from this gro file make an index file of protein without hydrogens
        - from the starting xtc of only protein make a xtc of protein without hydrogens
    """

    if os.path.exists(f"{parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}.tpr")==True:


        #convert tpr to create one with only protein
        try:
            p1 = subprocess.Popen(["echo", "1"],stdout=subprocess.PIPE)    #https://stackoverflow.com/questions/13332268/how-to-use-subprocess-command-with-pipes
            subprocess.run(shlex.split(parameters['gmx_path']+f" convert-tpr -s {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}.tpr -o {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_protein.tpr"), check = True, stdin=p1.stdout)
        except FileNotFoundError as exc:
            print_error(f"Process failed because the executable could not be found.\n{exc}")
        except subprocess.CalledProcessError as exc:
            print_error("convert-tpr doesn't work")
            sys.exit()


        #convert tpr to create one with only protein and no H
        try:
            p1 = subprocess.Popen(["echo", "2"],stdout=subprocess.PIPE)
            subprocess.run(shlex.split(parameters['gmx_path']+f" convert-tpr -s {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}.tpr -o {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.tpr"), check = True, stdin=p1.stdout)
            #subprocess.run(shlex.split(parameters['gmx_path']+" convert-tpr -s MD/"+parameters['md_new_name']+"/"+parameters['md_new_name']+".tpr -o MD/"+parameters['md_new_name']+"/"+parameters['md_new_name']+"_noH.tpr"), check = True)
        except FileNotFoundError as exc:
            print_error(f"Process failed because the executable could not be found.\n{exc}")
        except subprocess.CalledProcessError as exc:
            print_error("convert-tpr doesn't work")
            sys.exit()


        #create a gro file from the tpr with no H (needed later)
        try:
            subprocess.run(shlex.split(parameters['gmx_path']+f" editconf -f {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.tpr -o {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.gro"), check = True)
        except FileNotFoundError as exc:
            print_error(f"Process failed because the executable could not be found.\n{exc}")
        except subprocess.CalledProcessError as exc:
            print_error("editconf doesn't work")
            sys.exit()

        #make index file for protein with no H            
        try:
            entry="echo 2; echo q"        #https://stackoverflow.com/questions/17742789/running-multiple-bash-commands-with-subprocess
            p1 = subprocess.Popen(entry,stdout=subprocess.PIPE, shell=True)
            subprocess.run(shlex.split(parameters['gmx_path']+f" make_ndx -f {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.gro -o {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.ndx"), check = True, stdin=p1.stdout)
        except FileNotFoundError as exc:
            print_error(f"Process failed because the executable could not be found.\n{exc}")
        except subprocess.CalledProcessError as exc:
            print_error("make_ndx doesn't work")
            sys.exit()

        #make an xtc file of protein with no H
        try:
            p1 = subprocess.Popen(["echo", "2"],stdout=subprocess.PIPE)
            subprocess.run(shlex.split(parameters['gmx_path']+f" trjconv -f {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}.xtc -s {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_protein.tpr -o {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.xtc -pbc mol "), check = True, stdin=p1.stdout)
        except FileNotFoundError as exc:
            print_error(f"Process failed because the executable could not be found.\n{exc}")
        except subprocess.CalledProcessError as exc:
            print_error("trjconv doesn't work")
            sys.exit()

    else:
        print_error(f"File {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}.tpr not found")
        sys.exit()

def Make_RC(parameters):

    """if flag --rc is selected, Luigi will call multi-eGO to make a random coil .top and .gro of the reference structure, which will be inserted in outputs directory"""

    if parameters['rc']=='rc':
        try:
            subprocess.run( shlex.split(f"python multiego.py --protein={parameters['protein_name']} --egos=rc")   ,check=True)
            print_message('Luigi succesfully made random coil')
        except FileNotFoundError as exc:
            print_error(f"Process failed because the executable could not be found.\n{exc}")
        except subprocess.CalledProcessError as exc:
            print_error("multiego.py doesn't work")
            sys.exit()
    elif parameters['rc']=='none':
        print_message('Luigi will not do the Random Coil')

    else:
        print_error('Error with random coil parameters definition')
        sys.exit()


def Make_MDMAT(parameters):

    """
    This will use mdmat to calculate the contacts matrix of the md_ensemble trajectory selected as input in the --md flag of luigi.y.
    The output are dm.xpm, nat-all.ndx and mat.dat wich will be moved to the md_ensembles working directory
    """

    #calculate mdmat
    try:
        p1 = subprocess.Popen(["echo", "0"],stdout=subprocess.PIPE)    #https://stackoverflow.com/questions/13332268/how-to-use-subprocess-command-with-pipes
        subprocess.run(shlex.split(f"gmx_mpi mdmat -f {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.xtc -s {parameters['md_input_path']}/{parameters['md_new_name']}/{parameters['md_new_name']}_noH.tpr" ), check = True,  stdin=p1.stdout)
    except FileNotFoundError as exc:
        print_error(f"Process failed because the executable could not be found.\n{exc}")
    except subprocess.CalledProcessError as exc:
        print_error("mdmat doesn't work")
        sys.exit()
    
    #move mdmat files to md_ensamble working directory
    source='dm.xpm'
    destination=f"{parameters['md_input_path']}/{parameters['md_new_name']}/{source}"
    if os.path.exists(f'./{source}'):
        shutil.move(source, destination)
    else:
        print_error(f"{source} does not exists or is not accessible")

    source='nat-all.ndx'
    destination=f"{parameters['md_input_path']}/{parameters['md_new_name']}/{source}"
    if os.path.exists(f'./{source}'):
        shutil.move(source, destination)
    else:
        print_error(f"{source} does not exists or is not accessible")

    source='mat.dat'
    destination=f"{parameters['md_input_path']}/{parameters['md_new_name']}/{source}"
    if os.path.exists(f'./{source}'):
        shutil.move(source, destination)
    else:
        print_error(f"{source} does not exists or is not accessible")



##############################################################################################################################################################################################################################################

#MAIN

##############################################################################################################################################################################################################################################






def main(argv):

    print_hello()

    parameters={
        'protein_file':'none',
        'protein_name':'none',
        'protein_path':'none',
        'gmx_path':'none',
        'remove': [],
        'to_remove_list':'none',
        'md_path':'none',
        'md_name':'none',
        'md_new_name':'none',
        'rc':'none',

        'ref_path':'none',
        'md_input_path':'none',
        'protein_input_path':'none',
        'input_name':'none'
    }
    
    try:
        opts, args = getopt.getopt(argv,"g:p:r:m:n:i:ch",
                                ["gmx=",  "protein=", "remove=", 'md=', 'md_name=','input_name=', 'rc' ,"help"])
    except getopt.GetoptError:
        print_error('luigi.py --protein=<protein.pdb or PATH/proetin.pdb>, --gmx ,--remove, --md, --md_name, --rc, --input_name')
        sys.exit(2)
    if(len(opts)==0):
        print_message('luigi.py --protein=<protein.pdb or PATH/proetin.pdb>, --gmx, --remove, --md, --md_name, --rc, --input_name')
        sys.exit()


    #Read options
    parameters=read_options(opts, parameters)
    
    #REFERENCE preparation
    if parameters['protein_name']!='none':

        #Make REFERENCE directory
        Make_REF_directory(parameters)

        #Copy_protein in REFERENCE directory
        Copy_protein(parameters)

        #Remove from name.pdb all unnecessary molecules and rename it name_clean.pdb
        Remove_molecules(parameters)

        #Use pdb2gmx to create .gro .top
        Make_topology(parameters)

        #make random coil
        Make_RC(parameters)

    else:
        print_message('No protein name was inserted. Luigi will not prepare REFERENCE and use pdb2gmx')

    #MD preparation
    if parameters['md_name']!='none':

        #Make MD directory
        Make_MD_directory(parameters)

        #Copy md files in the correct directory
        Copy_md(parameters)

        #Make xtc file of only the protein without the hidrogens
        Make_xtc_noH(parameters)

        #Make the MDMAT of the md_ensemble selected as input in the --md flag of luigi.y
        Make_MDMAT(parameters)
    


    else:
        print_message('Luigi will not prepare MD')

    print_goodbye()
                

"""TO DO:  fare un clean degli input e degli output. fare il manuale. """
"""testare il copiaggio del force field """
"""trovare il modo di usare bene il comando gmx_mpi per mdmat"""
                



if __name__ == "__main__":
   main(sys.argv[1:])
