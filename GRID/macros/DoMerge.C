#include <TGrid.h>
#include <TFileMerger.h>

void DoMerge(
  const char* output="LHC15f_NewEvSel.root", 
  const char* path="/alice/cern.ch/user/f/fercoles/NOMECARTELLA/OutputDir/*"
  ) {

  system(Form("alien_ls -b %s/AnalysisResults.root |awk \'{print \"alien://\"$2}\' >listtobemerged",path));

  TGrid::Connect("alien://");
  FILE *f = fopen("listtobemerged","r");

  TFileMerger m(kFALSE);
  m.OutputFile(output);

  Int_t i=0;
  char nome[1000];
  while (fscanf(f,"%s",nome)==1) {
    m.AddFile(nome);
    i++;
  }
  if (i)
    m.Merge();
}

