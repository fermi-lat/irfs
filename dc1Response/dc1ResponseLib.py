#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['dc1Response'])
    env.Tool('st_facilitiesLib')
    env.Tool('irfInterfaceLib')
    env.Tool('irfUtilLib')
    env.Tool('astroLib')
    env.Tool('tipLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['f2cLibs'])

def exists(env):
    return 1
