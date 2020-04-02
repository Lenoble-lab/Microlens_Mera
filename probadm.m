% probabilite de masse des lentilles.

function pm = probadm(m)

global normfmdm
 
pm=fmdm(m)./normfmdm;