for i in range(1000):
    string="python add_systematics_13TeV_run2.py "+str(i)
    if i%7!=6:
      string+=" &"
    print string
