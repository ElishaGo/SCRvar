import subprocess
import sys
import os

file_dir = os.path.dirname(__file__)


def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', '-C', file_dir, 'rev-parse', 'HEAD']).decode('ascii').strip()


def get_git_branch_name() -> str:

    subprocess_res =  subprocess.check_output(['git', '-C', file_dir, 'branch', '--show-current'])
    branch_name = subprocess_res.decode('ascii').strip()
    return branch_name


def save_how_to(out_dir: str, sub_cmd_str: str = ''):
    """
    Saves the cml with git revision and in a txt file in the out_dir
    :param out_dir:
    :param sub_cmd_str: sub-command to append to 'how_to', e.g. how_to_eval.txt, how_to_train.txt (clipped at 20 chars)
    :return:
    """
    if len(sub_cmd_str) > 20:
        sub_cmd_str = sub_cmd_str[:20]
    with open(os.path.join(out_dir, 'how_to'+sub_cmd_str+'.txt'), 'wt') as f:
        f.write(' '.join(['python'] + sys.argv + ['\n']))
        f.write('git revision: ' + get_git_revision_hash() + '\n')
        f.write('git branch at time of execution: ' + get_git_branch_name() + '\n')
