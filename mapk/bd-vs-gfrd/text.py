#from paraview import servermanager
from paraview import simple


text_source = servermanager.sources.Text()
servermanager.Register(text_source, registrationName='text_source')
rep = simple.Show(text_source)

#rep = servermanager.CreateRepresentation(text_source, rv)
#servermanager.ProxyManager().RegisterProxy("representations",
#"my_representation%d" % 1, rep)

text_script = servermanager.sources.ProgrammableSource()
text_script.add_attribute('connection', servermanager.ActiveConnection)
text_script.Script = """
try:
    self.counter
except:
    self.counter = 0

if servermanager.ActiveConnection:
    print 'connection exists'
    #print servermanager.ActiveConnection
else:
    print 'create connection'
    #servermanager.Connect()

try:
    print self.connection
except:
    print 'no self.connection'
    pass

'''
rv = None
#rvs = servermanager.GetRenderViews()
#if len(rvs) > 0:
#  rv = rvs[0]
#else:
rv = servermanager.CreateRenderView()
'''

pxm = servermanager.ProxyManager()
text_source = pxm.GetProxy('sources', 'text_source')

text_source.Text = 'Ok' + str(self.counter)
self.counter += 1

'''
if rv.ViewTime.GetData() < 0.5:
    text_source.Text = '< 0.5'
    print '< 0.5'
else:
    text_source.Text = '> 0.5'
    print '> 0.5'
'''

"""
text_script.UpdatePipeline()

servermanager.Register(text_script, registrationName='text_script')


