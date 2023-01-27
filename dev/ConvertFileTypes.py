## Code to manipulate files (LGF)

## Grava preços e cashflow em excel a partir dos aruivos npy originais

gravaPlanilha = True

if gravaPlanilha:
    for date in example_dates:
        B=np.load(dir_data+'price_{}.npy'.format(date))
        print(B[:5])
        C=sps.load_npz(dir_data+'cashflow_{}.npz'.format(date)).toarray()
        print(C[:5])
        dados = np.load(dir_data+'cashflow_{}.npz'.format(date), allow_pickle=True)
                
        M=B.shape[0]
        print(M)
        
        print('Date: {}; Number of securities: {}'.format(date, M))
    
        Bdf = pd.DataFrame(B)
        Bdf.to_excel(dir_data+'preco_{}.xlsx'.format(date), sheet_name=date, index=False, header=False)
    
        Cdf = pd.DataFrame(C)
        Cdf.to_excel(dir_data+'fluxo_{}.xlsx'.format(date), sheet_name=date, index=False, header=False)
    
   

## Converte arquivos Excel para npz/npy

lePlanilha = True

if lePlanilha:
    for date in example_dates:
        Bdf = pd.read_excel(dir_data+'preco_{}.xlsx'.format(date), header=None, engine="openpyxl") 
        print(Bdf.head())
        B = np.array(Bdf[0])
        print(B.shape)
        np.save(dir_data+'preco_{}.npy'.format(date), B)
        
        Cdf = pd.read_excel(dir_data+'fluxo_{}.xlsx'.format(date), header=None, engine="openpyxl") 
        print(Cdf.head())
        C = np.asarray(Cdf)
        sps.save_npz(dir_data+'fluxo_{}.npz'.format(date), sps.csr_matrix(C))
        
        M=B.shape[0]
        print(M)

