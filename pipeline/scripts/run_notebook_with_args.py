import json
import subprocess
import sys

def run_notebook(notebook_file, output_file, *args):
    '''Pass arguments to a Jupyter notebook, run it and convert to html.
    From https://felixchenier.uqam.ca/a-better-way-to-run-a-jupyter-notebook-with-arguments/
    '''
    # Create the arguments file
    with open('arguments.json', 'w') as f:
        json.dump(args, f)
    # Run the notebook
    subprocess.call([
        'jupyter-nbconvert',
        '--execute',
        '--to', 'html',
        '--output', output_file,
        notebook_file])
    
def main():
    # First argument is the notebook path, second is output file path
    notebook_file = sys.argv[1]
    output_file = sys.argv[2]
    run_notebook(notebook_file, output_file, *sys.argv[3:])

if __name__ == "__main__":
    main()