bg_color white
set ray_opaque_background, off
set orthoscopic, 1
set antialias, 2
set ray_shadows, 0
set ray_trace_fog, 0
set cartoon_fancy_helices, 1
set cartoon_fancy_sheets, 1
set stick_radius, 0.10
set sphere_scale, 0.30
set label_size, 9
set label_font_id, 7
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.3

load ~/projects/rtc-pan-coronavirus/03-virtual-screening/NSP12-NSP7_3/receptor_NSP12-NSP7_3.pdb, NSP12_NSP7

hide everything, NSP12_NSP7
show cartoon, NSP12_NSP7
color gray85, NSP12_NSP7
set cartoon_transparency, 0.65, NSP12_NSP7

select nsp7, NSP12_NSP7 and chain C
color forest, nsp7
set cartoon_transparency, 0.10, nsp7

select anchor, NSP12_NSP7 and chain A and resi 440
show sticks, anchor
show spheres, anchor and name CA
set sphere_scale, 0.55, anchor and name CA
color red, anchor
label anchor and name CA, "PHE440"

select cluster, NSP12_NSP7 and chain A and resi 412+415+420+442
show sticks, cluster
color tv_orange, cluster
label cluster and name CA, "%s%s" % (resn, resi)

select triad, NSP12_NSP7 and chain A and resi 618+759+760
show sticks, triad
show spheres, triad and name CA
set sphere_scale, 0.40, triad and name CA
color marine, triad
label triad and name CA, "%s%s" % (resn, resi)

distance d_ASP618, anchor and name CA, (NSP12_NSP7 and chain A and resi 618 and name CA)
distance d_SER759, anchor and name CA, (NSP12_NSP7 and chain A and resi 759 and name CA)
distance d_ASP760, anchor and name CA, (NSP12_NSP7 and chain A and resi 760 and name CA)
distance d_NSP7,   anchor and name CA, (NSP12_NSP7 and chain C and resi 4  and name CA)

set dash_color, gray60, d_ASP618
set dash_color, gray60, d_SER759
set dash_color, gray60, d_ASP760
set dash_color, forest, d_NSP7
set dash_width, 2.5
set dash_gap,   0.25
set dash_radius, 0.035
set label_size, 8

mkdir ~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/pymol-sessions
mkdir ~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/visuals

orient NSP12_NSP7
zoom NSP12_NSP7, 5
save ~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/pymol-sessions/step01_spatial_map.pse
ray 2400, 1800
png ~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/visuals/Fig_Step01_3D_overview.png, dpi=300

zoom anchor, 18
turn y, 25
turn x, -10
ray 2400, 1800
png ~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/visuals/Fig_Step01_3D_anchor_zoom.png, dpi=300

print("Done")
