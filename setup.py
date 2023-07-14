from distutils.core import setup
setup(
  name = 'ToyMaker',         
  packages = ['ToyMaker'],
  version = '0.2.5',  
  license='MIT',    
  description = 'This library is designed to facilitate research of stochasticity effects in genetic circuits',
  author = 'Sergio Andrés Pachón Dotor',
  author_email = 'sap9827@gmail.com',
  url = 'https://github.com/SergioPachonDotor',
  download_url = 'https://github.com/SergioPachonDotor/ToyMaker/archive/refs/tags/v_0.2.5.tar.gz',
  keywords = ['Stochastic process', 'Genetic circuits', 'Simulation', 'Synthetic Biology', 'Systems Biology'],
  install_requires=[     
          'numpy',
          'numba',
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # "3 - Alpha", "4 - Beta" or "5 - Production/Stable"
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8', 
    'Programming Language :: Python :: 3.9', 
    'Programming Language :: Python :: 3.10',
  ],
)
