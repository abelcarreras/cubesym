try:
    from setuptools import setup, Extension
    use_setuptools = True
    print("setuptools is used.")
except ImportError:
    from distutils.core import setup, Extension
    use_setuptools = False
    print("distutils is used.")

setup(name='cubesym',
      version=1.0,
      description='cubesym module',
      author='Abel Carreras',
      url='https://github.com/abelcarreras/cubesym',
      author_email='abelcarreras83@gmail.com',
      packages=['cubesymapi'],
      scripts=['scripts/cubesym',
               'scripts/create_shapes',
               'scripts/density_example'],
      requires=['scipy', 'numpy', 'matplotlib'],
      license='MIT License'
      )
