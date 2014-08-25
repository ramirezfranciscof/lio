# 0.002394810

runprog=../lio/liosolo/liosolo
 echo 'estas son las bibliotecas que voy a usar'
 echo 'CHEQUEARLAS!!!!!!!!!!'
ldd ${runprog}
${runprog} -i olestra.in -b 631  -c oles.xyz -v > salida

