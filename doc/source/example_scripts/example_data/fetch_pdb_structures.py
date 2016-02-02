targets = ["7ODC","2QUD","3T47","4C1A","3ZSJ","4I4T","3H6R","3BSO",
           "3M0Z","2GPI","3MHP","3V19","2IXS","4G6T","3ZG9","3WSG",
           "1JH6","2XZE","2OLR","1TIF","3C4S","2CVI","4EET","1LM5",
           "1MJ5","2FJR","2D3G","3ZNV","3WA2","3WU2","4M7R","2PR7",
           "3FTJ","1KD8","1HBN","4TRK","2GB4","3HNY","4R7Q","1EAQ",
           "2O1Q","4DX5","1XAK","5CSM","2XWP","2UWA","3SQZ","4MT8",
           "3HTU","1CRN","3EML"]


for t in targets:
  prot = io.LoadPDB(t,remote=True)
  io.SavePDB(prot, t+".pdb")

