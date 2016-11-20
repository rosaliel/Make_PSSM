import re
from pprint import pprint

PATTERN = re.compile((
     '^ATOM|^HETATM.*ALA|^HETATM.*ARG|^HETATM.*ASN|^HETATM.*ASP|^HETATM.*CSH|'
     '^HETATM.*CYS|^HETATM.*GLN|^HETATM.*GLU|^HETATM.*GLY|^HETATM.*HIS|'
     '^HETATM.*ILE|^HETATM.*LEU|^HETATM.*LYS|^HETATM.*MET|^HETATM.*MSE|'
     '^HETATM.*ORN|^HETATM.*PHE|^HETATM.*PRO|^HETATM.*SER|^HETATM.*THR|'
     '^HETATM.*TRP|^HETATM.*TYR|^HETATM.*VAL|^HETATM.*ACE|^HETATM.*FOR|'
     '^HETATM.*ABA|^HETATM.*BOC|^HETATM.*BMT|^HETATM.*SAR|^HETATM.*MLE|'
     '^HETATM.*MVA|^HETATM.*IVA|^HETATM.*DFO|^HETATM.*NME|^HETATM.*AHT|'
     '^HETATM.*PTR|^HETATM.*PCA|^HETATM.*HYP|^HETATM.*INI|^HETATM.*NLE|'
     '^HETATM.*TYS|^HETATM.*CGU|^HETATM.*STA|^HETATM.*ILG|^HETATM.*OCS|'
     '^HETATM.*KCX|^HETATM.*SAH|^HETATM.*SAM|^HETATM.*SEP|^HETATM.*LLP|'
     '^HETATM.*5HP|^HETATM.*CSO|^HETATM.*ETA|^HETATM.*TFA|^HETATM.*ANI|'
     '^HETATM.*MPR|^HETATM.*DAM|^HETATM.*ACB|^HETATM.*ADD|^HETATM.*CXM|'
     '^HETATM.*DIP|^HETATM.*BAL|AL  |^HETATM.*AS |^HETATM.*BA|^HETATM.*BR|'
     '^HETATM.*CA  |^HETATM.*CD|^HETATM.*CL|^HETATM.*CO |^HETATM.*CS|'
     '^HETATM.*CU|^HETATM.*FE|^HETATM.*HG|^HETATM.*KR|^HETATM.*LA|^HETATM.*LI|'
     '^HETATM.*MG|^HETATM.*MN|^HETATM.*NA|^HETATM.*NI|^HETATM.*PB|^HETATM.*PR|'
     '^HETATM.*PT|^HETATM.*SE|^HETATM.*SM|^HETATM.*ZN'
     ))

def matching_lines(lines, pattern, invert=False):
    """Returns lines matching the pattern, or not matching in case 
    invert=True.
    """
    if invert:
        return [line for line in lines if not re.search(pattern, line)]
    else:
        return [line for line in lines if re.search(pattern, line)]

def sub(lines, pattern, repl):
    """"""
    return [re.sub(pattern, repl, line) for line in lines]
    
def clean_pdb(source, output):
    """"""
    # choose lines to correct
    lines = open(source, 'r').readlines()
    lines = matching_lines(lines, PATTERN) # grep all lines to consider later
    lines = matching_lines(lines, '[B-D][A-Z][A-Z][A-Z]', invert=True)
    lines = matching_lines(lines,
                           'SEQRES|REMARK|HELIX|SHEET|TITLE|JRNL|SITE|LINK',
                           invert=True)
    lines = matching_lines(lines, 
                           'HETATM.*CX.*KCX|HETATM.*OQ1.*KCX|HETATM.*OQ2.*KCX',
                           invert=True)
    
    # correct lines with weird AAs
    lines = sub(lines, '(HETATM)(.*)(KCX)', 'ATOM  \\2LYS')
    lines = sub(lines, '(HETATM)(.*)(MSE)', 'ATOM  \\2MET')
    lines = sub(lines, '(ATOM   .*)SE (  MET)', '\\1 SD\\2')
    lines = sub(lines, '(HETATM)(.*)(CSO)', 'ATOM  \\2CYS')
    lines = sub(lines, '(HETATM\(.*)(CSH)', 'ATOM  \\2CYS')
    
    # remove any non atom lines
    lines = matching_lines(lines, '^ATOM')   
    
    # remove Alternate location indicator
    lines = [line[:16] + ' ' + line[17:77] + '\n' for line in lines] 
    
    open(output, 'w').write(''.join(lines))

    
    
if __name__ == '__main__':
    clean_pdb('5t5i.pdb','m')
    