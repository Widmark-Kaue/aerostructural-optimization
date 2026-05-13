'''
Module with Tapenade interface

Ney Rafael Secco
2021-04-04
'''

import os
import fileinput
import argparse

#===================================================
# USER INPUTS
'''
This function takes files from the current directory and submits them to
the local Tapenade interface for automatic differentiation.

INPUTS (with examples)

files: list of strings -> List of filenames that should be differentiated.
files = ['a_module.f90',
         'b_module.f90',
         'c_module.f90']

top_routine: string -> Name of the top routine that will be differentiated.
top_routine = 'main'

input_vars: string -> String with input variable names separated by whitespaces.
input_vars = 'X Y Z'

output_vars: string -> String with output variable names separated by whitespaces.
output_vars = 'f g h'
'''

# List of files that will be differentiated
files = ['???']

# Name of top routine
top_routine = '???'

# Output variables (separate with whitespace)
output_vars = '???'

# Inputs variables (separate with whitespace)
input_vars = '???'

#===================================================

# GET DIFF MODE FROM COMMAND LINE ARGUMENT

# Create parser
parser = argparse.ArgumentParser(description='Send files to the Tapenade web interface for automatic differentiation.')
parser.add_argument('mode', action='store',
                    help='differentiation mode: f (forward)) or r (reverse)')
args = parser.parse_args()

if args.mode == 'f':
    mode = 'forward'
elif args.mode == 'r':
    mode = 'reverse'
else:
    raise ValueError('Differentiation mode not recognized. Use either f or r.')

#===================================================

def diff(files,
         top_routine,
         input_vars,
         output_vars,
         mode,
         output_directory = 'TapenadeResults',
         real4toreal8 = True,
         diff2mode = True):

    '''
    This function takes files from the current directory and submits them to
    the local Tapenade interface for automatic differentiation.
    
    INPUTS (with examples)
    
    files: list of strings -> List of filenames that should be differentiated.
    files = ['a_module.f90',
             'b_module.f90',
             'c_module.f90']
    
    top_routine: string -> Name of the top routine that will be differentiated.
    top_routine = 'main'
    
    input_vars: string -> String with input variable names separated by whitespaces.
    input_vars = 'X Y Z'
    
    output_vars: string -> String with output variable names separated by whitespaces.
    output_vars = 'f g h'
    
    mode: string -> Either 'forward' or 'reverse'.
    mode = 'forward'
    
    output_directory: string -> Directory where differentiated codes
                                will be stored.
    
    real4toreal8: boolean -> Flag to replace the PUSH/POP functions inserted
                         in the reverse code from REAL4 to REAL8.
    
    diff2mode: boolean -> Flag to replace '_diff' in the differentiated
                          module names by either '_d' or '_b'.
    '''

    # Generate directory name to extract differentiated files
    if mode == 'forward':
        dirname = output_directory + '_d'
    elif mode == 'reverse':
        dirname = output_directory + '_b'

    # Check if we already created the folder with differentiated files
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # BUILD THE TAPENADE COMMAND

    # Define routines and variables that will be differentiated
    diff_comm = '-head "%s(%s)\\(%s)"'%(top_routine,input_vars,output_vars)

    # Differentiation mode flag
    if mode == 'forward':
        mode_flag = '-d'
    elif mode == 'reverse':
        mode_flag = '-b'

    # Directory where files will be stored
    output_flag = '-O '+dirname

    # List of files that will be differentiated
    files_comm = ' '.join(files)

    # Run the command line
    os.system('tapenade %s %s %s %s'%(diff_comm,
                                      mode_flag,
                                      output_flag,
                                      files_comm))

    # Check if we need to replace real4 by real8 in reverse differentiated files
    if (mode == 'reverse') and real4toreal8:

        # Print log
        print('Replacing REAL4 by REAL8')

        # Get list of files in the extracted directory
        files = os.listdir(dirname)

        # Loop over the differentiated file names
        for ff in files:

            # Check if this is a differentiated file
            if '_b.f90' in ff:

                # Get full filename
                ff_full = os.path.join(dirname,ff)

                # Replace real4 by real8
                for line in fileinput.input([ff_full], inplace=True):
                    print(line.replace('REAL4', 'REAL8'), end='')

    # Check if we need to replace diff by the differention mode in module names
    if diff2mode:

        # Print log
        print('Replacing DIFF in module names')

        # Get list of files in the extracted directory
        files = os.listdir(dirname)

        # Set file differentiation tag
        if mode == 'forward':
            filetag = '_d.f90'
            modetag = '_D'
        elif mode == 'reverse':
            filetag = '_b.f90'
            modetag = '_B'

        # Loop over the differentiated file names
        for ff in files:

            # Check if this is a differentiated file
            if filetag in ff:

                # Get full filename
                ff_full = os.path.join(dirname,ff)

                # Replace real4 by real8
                for line in fileinput.input([ff_full], inplace=True):
                    print(line.replace('_DIFF', modetag), end='')


#===================================================

# Now that the function is defined, we can execute it
# with the user's parameters
diff(files,
     top_routine,
     input_vars,
     output_vars,
     mode)