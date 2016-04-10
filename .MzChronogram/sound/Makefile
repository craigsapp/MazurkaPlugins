
sweep:

	/project/websites/sig/bin/osc -d 5  \
	-f "0 43.06640625 1 43.06640625 2 215.33203125 3 215.33203125 4 430.6640625 5 430.6640625" \
        -a "0 0 10 0.5 90 0.5 100 0" sweep43-430.wav


analysis:

	../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mztimeogram \
	sweep43-430.wav 0 | sed 's/.*: //' > sweep43-430.txt

clean:

	rm -i sweep43-430.wav sweep43-430.txt


chromo:

	smf -CT0 -o chromo.mid chromo.krn

	timidity -c /usr/share/timidity/eawpatches/timidity.cfg  -OwM \
	 	-o chromo.wav chromo.mid


chromoanalysis:

	../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mztimeogram \
	chromo.wav 0 | sed 's/.*: //' > chromo.txt

eightfreq:

	/project/websites/sig/bin/osc -d 9 \
	-f "0 430.6640625 100 430.6640625 101 436.04736328125 200 436.04736328125 201 441.4306640625 300 441.4306640625 301 446.81396484375 400 446.81396484375 401 452.197265625 500 452.197265625 501 457.58056640625 600 457.58056640625 601 462.9638671875 700 462.9638671875 701 468.34716796875 800 468.34716796875 801 473.73046875 900 473.73046875" \
	-a "0 0 3 0.24 97 0.24 100 0" eightfreq.wav


eightfreqanalysis:

	../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mztimeogram \
	eightfreq.wav 0 | sed 's/.*: //' > eightfreq.txt
