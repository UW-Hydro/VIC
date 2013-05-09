BEGIN {
    dp=42
    yp=30
    mp=34
    pp=22
    ip=6
    lmo=0
    monc=0
    dys="31:29:31:30:31:30:31:31:30:31:30:31"
    split(dys,ndy,":")
    fact=1.0
    rec1=0
    }
    {
    param=substr($0,pp,4)
    if(NR>2 && param=="PRCP"){
	id=substr($0,ip,6)
	if(id!=999999){
		rec1++
		if(param=="PRCP") fact=100.0
		yr=substr($0,yp,4)
		mon=substr($0,mp,2)
		yr+=0
		mon+=0
                idstr=id
		id+=0
		if(rec1==1) {
		    yc=yr
        	    idc=id
		    printf("%6s %4s\n",idstr,yc)
        	}
                if(id != idc || yr != yc || mon != monc) {
		 if(id != idc || yr != yc ) {
		    if(lmo > 0) {
			for(j=lmo+1;j<=12;j++){
			    for(i=1;i<=ndy[j];i++) printf("-99 ")
			}
		    }
		    printf("\n")
		    idc=id
		    yc=yr
		    lmo=0
		    printf("%6s %4s\n",idstr,yc)
		 }
		 if(mon > lmo+1) {
		    for(j=lmo+1;j<mon;j++) {
			for(i=1;i<=ndy[j];i++) printf("-99 ")
	            }
		    lmo=mon-1
		 }
		 if(mon==lmo+1){
		    for(i=1;i<=ndy[mon];i++) {
			dt=substr($0,dp+(16*(i-1)),6)
			dt*=1.0
			if(dt<-99998 || dt>99998) printf("-99 ")
			else printf("%.2f ",dt/fact)
		    }
		 }
		if(mon>lmo) {lmo=mon}
		if(lmo==12)lmo=0
                monc=mon
                }
    	}
    }
    }
    END{
    if( mon < 12 ) {
	for(j=mon+1;j<=12;j++){
		    for(i=1;i<=ndy[j];i++) printf("-99 ")
	}
    }
    printf("\n")
    }
