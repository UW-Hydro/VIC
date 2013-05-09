      program argtest

      character tpsw*2
      character infofile*80
      character maskfile*80
      character inputfile*80
      character outputfile*80


      call getarg(1,tpsw)
      call getarg(2,infofile)
      call getarg(3,maskfile)
      call getarg(4,inputfile)
      call getarg(5,outputfile)

      if(tpsw.eq.''.or.infofile.eq.''.or.maskfile.eq.''.or.inputfile.
     $eq.''.or.outputfile.eq.'') then
      print*,'Usage:'
      print*,'regrid [t/p] [infofile] [maskfile] [input] [output]'
      stop
      end if

      print*, tpsw,infofile,maskfile,inputfile,outputfile

      stop
      end
      
