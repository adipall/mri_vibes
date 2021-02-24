#module load sierra
#source /data1/opt/sierra/install/sierra_init.sh
#source /etc/profile

#updated 10/21/2020
source /opt/intel/compilers\_and\_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
export SIERRA_INSTALL_DIR=/data1/opt/sierra.4.52.2/install
#source \SIERRA\_INSTALL\_DIR/sierra\_init.sh
echo $SIERRA_INSTALL_DIR
source /data1/opt/sierra.4.52.2/install/sierra\_init.sh
source /etc/profile
