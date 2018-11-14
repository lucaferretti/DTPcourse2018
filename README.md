# Start of the session:

first open a terminal, then connect to the Pirbright cluster via SSH

> ssh DTCstudents@mallorn.pirbright.ac.uk

and use the password that was given during the course. Then ask one of the demonstrators for the numeric code to enter.


The first time that you connect, a few preparatory steps to work remotely:

choose a (unique!) username, containing only letters and numbers and underscores, then type

> tmux new -s MYUSERNAME

where MYUSERNAME is the username you chose. This will open a virtual session on the cluster that will avoid losing your work if you lose your connection!

Then create your own directory

> mkdir MYUSERNAME

copy all data there:

> cp Data/* MYUSERNAME/

then move there to work

> cd MYUSERNAME

Please do not work outside this directory!

Finally, reserve a computational unit for yourself by running

> srun -c 1 --tasks-per-node 1 --pty bash -i

This will open you a shell on a computational node.

From now on, please avoid typing "exit" or anything similar. In case you want to get out, just close the window of the terminal.


After the first time, you can simply open a terminal again, reconnect via SSH as described above, then type

> tmux new -t MYUSERNAME

and you'll get to the same point you were when you left.

----------------------


