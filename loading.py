import pymol

x_length = 1500
y_length = 800

obj_protein = '3vi8'
obj_pah = 'pah_aligned'

color_protein = 'marine'
color_ligand = 'red'
color_pah = 'limon'
sel_protein = 'prot'
sel_ligand = 'lig'
sel_pah = 'pah_aligned'
sel_backbone = 'n. C+CA+O+N'

attribute_cartoon = 'cartoon'
attribute_stick = 'stick'
attribute_surf = 'surface'

whole_view = '''\
      0.416837215,   -0.585738242,   -0.695094168,\
      0.402514488,   -0.566689312,    0.718919933,\
     -0.814999342,   -0.579462051,   -0.000447881,\
      0.000000000,    0.000000000, -132.898391724,\
      5.669492722,    1.393989563,    0.986843109,\
    -172.826446533,  438.675079346,  -20.000000000'''

pocket_view_1 = '''\
      0.570410669,   -0.507711053,    0.645637393,\
     -0.150874689,    0.707922995,    0.689984441,\
     -0.807378590,   -0.490983635,    0.327207208,\
     -0.000129877,    0.000011414,  -68.655761719,\
      9.738702774,    4.544628620,   -6.171709061,\
    -12.349261284,  152.696075439,  -20.000000000'''

cutaway_view_1 = '''\
     0.260635376,    0.221995786,   -0.939556897,\
     0.746257722,   -0.663764596,    0.050179087,\
    -0.612505674,   -0.714234769,   -0.338663787,\
    -0.000129877,    0.000011414,  -68.655761719,\
     9.738702774,    4.544628620,   -6.171709061,\
    66.461433411,   73.885368347,  -20.000000000'''

cutaway_view_1_1 = '''\
     0.422939479,    0.185698941,   -0.886915028,\
     0.713203490,   -0.672000647,    0.199398950,\
    -0.558979928,   -0.716889381,   -0.416653514,\
    -0.000074793,   -0.000075666,  -69.355690002,\
    10.387354851,    5.356313705,   -5.510336399,\
    65.561447144,   73.885368347,  -20.000000000'''

cutaway_view_2 = '''\
     0.969169736,    0.018895648,   -0.245457664,\
    -0.110522844,   -0.857491910,   -0.502406895,\
    -0.219995558,    0.514051616,   -0.829003334,\
    -0.000949420,    0.000929267,  -58.990161896,\
     9.557539940,    6.241314411,   -7.982754707,\
    54.210052490,   69.668014526,  -20.000000000'''

residue_view = '''\
    -0.267488927,    0.833759785,    0.482891977,\
    -0.564587951,   -0.541739941,    0.622632027,\
     0.780755043,   -0.106111921,    0.615672112,\
     0.001483641,   -0.000500267,  -58.785842896,\
     9.956514359,    7.084040642,   -9.176513672,\
  -971.948364258, 1095.256713867,  -20.000000000'''

def setup():
    cmd.reinitialize()
    cmd.load('{}.pdb'.format(obj_protein))
    cmd.load('{}.pdb'.format(obj_pah))
    cmd.reset() # reset view
    hide_all()
    util.performance(0) # maximum quality
    cmd.bg_color('white') # white background
    cmd.set('cartoon_fancy_helices','1',quiet=0) # fancy helixes

    cmd.clip('atoms','200','*')
    cmd.select(sel_protein,'(br. {0}) & {1}'.format(sel_backbone,obj_protein))
    cmd.select(sel_ligand,'{0}  & !{1}'.format(obj_protein,sel_protein))

    cmd.set('surface_color_smoothing_threshold', 0.7)
    cmd.set('surface_color_smoothing',20)
    cmd.set('transparency_mode',3,quiet=0)
    cmd.set('surface_quality',2)
    cmd.set('surface_cavity_mode', 2)
    cmd.set('cavity_cull', 150)    

def hide_all():
    cmd.hide('everything', '*')

def show_pah():
    cmd.show(attribute_stick, sel_pah)
    cmd.hide('(hydro & ' + sel_pah + ')')
    cmd.color(color_pah, sel_pah)
    
def show_ribbon(transp=0):
    cmd.show(attribute_cartoon, sel_protein)
    cmd.color(color_protein, sel_protein)
    cmd.set('cartoon_transparency',transp)

def show_ligand():
    cmd.show(attribute_stick, sel_ligand)
    cmd.color(color_ligand, sel_ligand)

def show_surface(transp=0.7):
    cmd.show(attribute_surf, sel_protein)
    cmd.set('transparency', transp, sel_protein)
    cmd.color(color_protein, sel_protein)
    # util.cnc(sel_protein)

def show_nearby_residues():
    rad = 6
    sel_nearby_atom = '{0!s} within {1!s} of {2!s}'\
                      .format(sel_protein,rad,sel_pah)
    sel_nearby_res = '(br. ({0})) & !({1}) & !({2})'\
                     .format(sel_nearby_atom,'n. C+O+N',sel_pah)
    cmd.select('{0!s}_of_Lig_Res'.format(rad),sel_nearby_res)
    sel_nearby_res = '{0!s}_of_Lig_Res'.format(rad)
    cmd.show(attribute_stick,sel_nearby_res)
    cmd.color(color_protein, sel_protein)
    util.cnc(sel_nearby_res)
    # cmd.label(sel_nearby_res, 'resi')
    cmd.label('''n. CA+C1*+C1' & (br. 6_of_Lig_Res)''','''"%s-%s"%(resn,resi)''')
    cmd.set('label_size',-0.5)
    cmd.set('label_color','black')
    cmd.set('label_outline_color','black')
    set_view(residue_view)
    # attribute_ribbon = 'ribbon'
    # cmd.show(attribute_ribbon, sel_nearby_res)
    
    
def save_image(name='UntiledImage',ray=1):
    print('Saving Image {0}_{1}.png ...'.format(obj_protein, name))
    if ray: cmd.ray()
    cmd.png(obj_protein + '_' + name)
    print('Saving Done.')

def print_sec_headder(number='?', title='Untitled Section'):
    print('#'*20 + ' Section ' + str(number)  + ': '  + title)
    
def set_view(view_str):
    cmd.set_view(view_str)
    cmd.deselect()
    
def main(section=0, render=0):
    sec_number = 1
    if not section == 0: sec_number = section

    if section == 0 or section == 1:
        print_sec_headder(sec_number,'Setup')
        setup()
        sec_number += 1

    if section == 0 or section == 2:
        print_sec_headder(sec_number,'Cartoon View')
        hide_all()
        show_ligand()
        show_ribbon()
        set_view(whole_view)
        if render: save_image('cartoon')
        sec_number += 1

    if section == 0 or section ==  3:
        print_sec_headder(sec_number,'Surface Pocket View')
        hide_all()
        show_ligand()
        show_surface()
        set_view(pocket_view_1)
        if render: save_image('surface')
        show_pah()
        if render: save_image('surface_PAH')
        sec_number += 1

    if section == 0 or section == 4:
        print_sec_headder(sec_number,'Pocket View With Cartoon')
        hide_all()
        show_ligand()
        show_ribbon(0.5)
        show_surface(0.3)
        set_view(whole_view)
        if render: save_image('cartoonAndPocket')
        show_pah()
        if render: save_image('cartoonAndPocket_PAH')
        sec_number += 1

    if section == 0 or section == 5:
        print_sec_headder(sec_number,'Cutaway Pocket View #1')
        hide_all()
        show_ligand()
        show_surface(0.3)
        set_view(cutaway_view_1_1)
        if render: save_image('cutaway1', 0)
        show_pah()
        if render: save_image('cutaway1_PAH', 0)
        sec_number += 1

    if section == 0 or section == 6:
        print_sec_headder(sec_number,'Cutaway Pocket View #2')
        hide_all()
        show_ligand()
        show_surface(0)
        set_view(cutaway_view_2)
        if render: save_image('cutaway2')
        show_pah()
        if render: save_image('cutaway2_PAH')
        sec_number += 1

    if section == 0 or section == 7:
        print_sec_headder(sec_number,'Nearby Residues')
        hide_all()
        show_ligand()
        show_nearby_residues()
        if render: save_image('residues')
        show_pah()
        if render: save_image('residues_PAH')
        sec_number += 1
                                                              
    print('*'*50+ ' '*2 + 'DONE' + ' '*2 + '*'*50)
            
    
cmd.extend('main', main)

