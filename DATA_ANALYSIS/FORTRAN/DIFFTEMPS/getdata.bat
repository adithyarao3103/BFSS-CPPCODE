FOR %%t IN (0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28) DO (
    FOR %%f IN (1 2 3 4) DO (
        scp draco:~/bfss_fortran/runs/difftemps32/%%t/%%t%%f.dat D:/DRACOFORTRAN/DIFFTEMPSNEW/
    )
)