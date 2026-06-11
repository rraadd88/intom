## creates virt. env. with two kernels one for each python and R 

env_name='intom'
## kernel names
kernel=$env_name
kernel_r=${env_name}_r

## package manager
pm='micromamba' # mamba/conda

## env 
com="
$pm create -n ${env_name} python=3.12 -y;
$pm activate ${env_name};
"
echo $com

## python packages
com="pip install .;
$pm install -c conda-forge rpy2 xz;
pip show rpy2;
"
echo $com

## py kernel
com="
pip install ipykernel;
python -m ipykernel install --user --name ${kernel} --display-name ${kernel};
"
echo $com

## R kernel
com="
$pm install -c r r-essentials r-base notebook r-irkernel r-repr r-irdisplay r-pbdzmq r-devtools r-biocmanager --channel-priority flexible -y;
$pm install bioconda::bioconductor-limma=3.62.1 -y;
R -e \"IRkernel::installspec(name='${kernel_r}',displayname='${kernel_r}')\";
"
echo $com

com="
sudo ln -s $(which R) /usr/local/bin/R
sudo ln -s $(which Rscript) /usr/local/bin/Rscript
"
echo $com