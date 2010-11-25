from paraview import servermanager
from paraview import simple

text = servermanager.sources.ProgrammableFilter()
filter.Initialize()
text.Script = """#from paraview import simple
'''
view = simple.GetActiveView()
if view.ViewTime < 0.5:
  print 'hoi'
else:
  print '> 0.5'
'''

"""
text.UpdatePipeline()

servermanager.Register(text, registrationName='text')


