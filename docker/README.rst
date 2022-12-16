# dockerfiles

.. code-block:: bash
		
    version=1.0b
    name=nanorms
    echo $name:$version
    
    docker build --pull -t lpryszcz/$name:$version .
    docker tag lpryszcz/$name:$version lpryszcz/$name:latest

    docker push lpryszcz/$name:$version && docker push lpryszcz/$name:latest

    # testing
    cd ~/test/modPhred/test;
    v=3.6.1; acc=PRJEB22772; docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/modphred-$v /opt/modPhred/run -f /data/ref/ECOLI.fa -o /data/modPhred.docker.$v/$acc -i /data/$acc/{MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145,MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711} -t4 --host 10.46.1.65 -p 5556

    v=5.0.11; acc=PRJEB22772; docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/modphred-$v /opt/modPhred/run -f /data/ref/ECOLI.fa -o /data/modPhred.docker.$v/$acc -i /data/$acc/{MARC_ZFscreens_R9.4_1D-Ecoli-run_FAF05145,MARC_ZFscreens_R9.4_2D-Ecoli-run_FAF05711} -t4 --host /usr/bin/guppy_basecall_server -c dna_r9.4.1_450bps_modbases_5mc_hac.cfg
