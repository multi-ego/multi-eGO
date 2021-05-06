awk -v atom=$i '{a1[NR]=$1; a2[NR]=$2} END {
      j=1;
      first=1;
      for(i=1;i<=NR;i++) {
        if(a1[i]==atom) {excl[j]=a2[i]; first=j; j++}
        else if(a2[i]==atom) {excl[j]=a1[i]; first=j; j++}
        else continue;
        for(k=1;k<=NR;k++) {
          second=1;
          if(a1[k]==excl[first]&&a2[k]!=atom) {excl[j]=a2[k]; second=j; j++}
          else if(a2[k]==excl[first]&&a1[k]!=atom) {excl[j]=a1[k]; second=j; j++}
          else continue;
          for(l=1;l<=NR;l++) {
            if(a1[l]==excl[second]&&a2[l]!=excl[first]) {excl[j]=a2[l]; j++}
            else if(a2[l]==excl[second]&&a1[l]!=excl[first]) {excl[j]=a1[l]; j++}
          }
        }
      }
      for(i=1;i<j;i++) print atom, excl[i]
     }' bond