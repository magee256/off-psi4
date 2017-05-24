
### Tcl script written to select conformers for some of diverse set to fulfill
#     specified test criteria, e.g, pi stacking.
#     Takes a mol2 file and outputs a filtered mol2 file with confs fitting the criteria.
#     The original oemol coordinates are retained.
#     
# Usage: vmdt -e file.tcl -args infile.mol2 outfile.mol2


set verbose 0
set lastF -1

### open data output file
set outmol2 [lindex $argv 1]
if {$verbose} {
    set outFile [open [lindex $argv 2] w]
    puts $outFile "#Conformer\t\tRing distance (A) "
}

mol new [lindex $argv 0] first 0 last $lastF waitfor all
set nframes [molinfo top get numframes]

for {set i 0} {$i < $nframes} {incr i} {

    ### Div_4 ================================================
    # use geometric center of rings (pi stacking)
    set vec1 [measure center [atomselect top "name C1 C2 N1 C5 N3 C6" frame $i]]
    set vec2 [measure center [atomselect top "name N2 C3 C4 N4 C7" frame $i]]
    set arodist [veclength [vecsub $vec1 $vec2]]
    # mark frames to delete, dist >= 5 angstroms
    if {$arodist >= 5} {lappend delist $i}
    # ========================================================

#    ### Div_5 ================================================
#    # use geometric center of rings (pi stacking)
#    set vec1 [measure center [atomselect top "name C1 C3 C4 C7 C8 C11" frame $i]]
#    set vec2 [measure center [atomselect top "name C2 C5 C6 C9 C10 C12" frame $i]]
#    set arodist [veclength [vecsub $vec1 $vec2]]
#    # mark frames to delete, dist >= 5 angstroms
#    if {$arodist >= 5} {lappend delist $i}
#    # ========================================================


    if {$verbose} {puts $outFile "$i\t$arodist"}

}

### loop thru flagged frames in reverse so as to not ruin count by deleting
foreach x [lreverse $delist] {
    animate delete beg $x end $x skip 0 0
}

if {$verbose} {close $outFile}
animate write mol2 $outmol2 beg 0 end -1 skip 1 0
exit
