import pytest


# Tests the generation of the symmetries from the oper_axpre of mmCIF and PDB files
@pytest.mark.parametrize('expression,idx_list', [
    ('(1-2)',[['1'],['2']]),
    ('(1-20)', [[str(x)] for x in range(1,21)]),
    ('1,3,5-15', [['1'],['3']]+[[str(x)] for x in range(5,16)]),
    ('(1-20)(21-25)(26,28)', [[str(i), str(j), str(k)] for i in range(1,21) for j in range(21,26) for k in [26, 28]])
  ])
def test_read_oper_expression(expression, idx_list):
    from supporting_functions import read_oper_expression
    mtx_ind = [str(x) for x in range(1,29)]
    exp_idx_list = read_oper_expression(expression, mtx_ind)
    assert exp_idx_list == idx_list

@pytest.mark.parametrize('elements_list,merged_list', [
    ([[0, 1], [1, 2], [2, 3], [1, 4], [5, 6]], [[0, 1, 2, 3, 4], [5, 6]]),
    ([[0], [1, 2, 3, 4, 5], [0, 1]], [[0, 1, 2, 3, 4, 5]]),
    ([[0, 2, 4], [1, 3, 5], []],  [[0, 2, 4], [1, 3, 5]])
  ])
def test_merge(elements_list, merged_list):
    from supporting_functions import merge
    exp_merged_list = merge(elements_list, sorted_lists=True)
    assert exp_merged_list == merged_list


@pytest.mark.parametrize('numbers_list,gcd', [
    ([1, 3, 5, 7], 1),
    ([2, 4, 6, 8], 2),
    ([18, 24, 48], 6),
    ([2,6,8,2,6,6,10,6,6,14,12], 2)
  ])
def test_gcd(numbers_list, gcd):
    from supporting_functions import find_gcd
    exp_gcd = find_gcd(numbers_list)
    assert exp_gcd == gcd


@pytest.mark.parametrize('pdbi', [
    ('3dh4'),
    ('1bcc'),
    ('1okc')
  ])
def test_mmCIF_format_parser(pdbi):
    import os
    import subprocess
    from supporting_functions import parse_mmCIF_format
    cif_filename = "{0}.cif".format(pdbi.upper())
    if not os.path.exists(cif_filename):
        subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(cif_filename)])
    n_loops = int(subprocess.run([f"grep 'loop_' {cif_filename} | wc -l"], shell=True, stdout=subprocess.PIPE).stdout.decode('utf8').split()[0])
    n_all_loops = int(subprocess.run([f"grep '^#\\s*$' {cif_filename} | wc -l"], shell=True, stdout=subprocess.PIPE).stdout.decode('utf8').split()[0])-1
    n_keys = int(subprocess.run([f"grep '^_' {cif_filename} | wc -l"], shell=True, stdout=subprocess.PIPE).stdout.decode('utf8').split()[0])
    all_loops, loop_keys = parse_mmCIF_format(cif_filename)
    assert len([x for x in all_loops if len(x)>2]) == n_loops
    assert len(all_loops) == n_all_loops
    assert len(loop_keys) == n_keys


@pytest.mark.parametrize('pdbi', [
    ('3dh4'),
    ('1bcc'),
    ('1okc')
  ])
def test_header_parsers(pdbi):
    import os
    import subprocess
    from supporting_functions import parse_mmCIF_format, parse_mmCIF_header, parse_PDB_header
    cif_filename = "{0}.cif".format(pdbi.upper())
    if not os.path.exists(cif_filename):
        subprocess.run(["wget", "https://files.rcsb.org/view/{0}".format(cif_filename)])
    all_loops, loop_keys = parse_mmCIF_format(cif_filename)
    mmCIF_header = parse_mmCIF_header(all_loops, loop_keys)
    operators = set([k for x in mmCIF_header[0]['biological_assemblies'] for y in x for z in y['operators'] for k in z])
    ktranslrots = set([x for x in mmCIF_header[0]['ktranslrots']])
    assert operators == ktranslrots



