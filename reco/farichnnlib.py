# -*- coding: utf-8 -*-
"""farichnnlib.py
"""

import sys, os, time, subprocess
import tensorflow as tf
import numpy as np
from math import pi
import pandas as pd
from sklearn.model_selection import train_test_split
import uproot3 as uproot
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['font.size'] = '16'
from particle import Particle
from iminuit import Minuit
from probfit import BinnedLH, gaussian, linear, AddPdfNorm
if 'useCuPy' in globals() and useCuPy:
  import cupy as cp


speedOfLight_mmperns = 299.792458 # мм/нс

def initcompstrategy(oncolab: bool):
  '''
  Инициализировать стратегию вычисления strategy в зависимости от доступной среды выполнения: CPU, GPU, TPU.
  '''
  strategy = None

  if oncolab:
    try:
      resolver = tf.distribute.cluster_resolver.TPUClusterResolver()
      tf.config.experimental_connect_to_cluster(resolver)
      # This is the TPU initialization code that has to be at the beginning.
      tf.tpu.experimental.initialize_tpu_system(resolver)
      strategy = tf.distribute.experimental.TPUStrategy(resolver)
      tf.compat.v1.disable_eager_execution()
    except Exception as e:
      print('Could not set TPU distribution strategy!')
      print(e)
      pass

  phys_devs = tf.config.experimental.list_physical_devices()
  phys_gpus = tf.config.experimental.list_physical_devices('GPU')
  phys_tpus = tf.config.experimental.list_physical_devices('TPU')

  if strategy is None and len(phys_gpus) > 0:
    strategy = tf.distribute.MirroredStrategy()
    try:
      if not oncolab:
        # Restrict TensorFlow to only allocate 1GB of memory on the first GPU if on local environment
        tf.config.experimental.set_virtual_device_configuration(phys_gpus[0],
                                                  [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=1024)])
      else:
        for gpu in phys_gpus:
          tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
      # Virtual devices must be set before GPUs have been initialized
      print(e)

  if strategy is None:
    strategy = tf.distribute.get_strategy()

  logical_devs = tf.config.experimental.list_logical_devices()
  logical_gpus = tf.config.experimental.list_logical_devices('GPU')
  logical_tpus = tf.config.experimental.list_logical_devices('TPU')

  for devs in ((phys_devs, 'devices'), (phys_gpus, 'GPUs'), (phys_tpus, 'TPUs')):
    print(f'{len(devs[0])} physical {devs[1]}:',' '.join([d.name for d in devs[0]]))
  for devs in ((logical_devs, 'devices'), (logical_gpus, 'GPUs'), (logical_tpus, 'TPUs')):
    print(f'{len(devs[0])} logical {devs[1]}:',' '.join([d.name for d in devs[0]]))

  return strategy

def show_uproot_tree(obj, maxkeylen=12, sep='/', indent=0):
  '''
  Отображение имеющихся ключей в ROOT-файле виде иерархии.
  '''
  width = maxkeylen+len(sep)
  startline = False
  if isinstance(obj, uproot.rootio.ROOTDirectory):
    print('TFile: '+obj.name.decode('utf-8'))
    startline = True
    indent = 2
  elif issubclass(type(obj), uproot.tree.TTreeMethods):
    print('TTree: '+obj.name.decode('utf-8'))
    startline = True
    indent = 4
  else:
    if len(obj.keys()) > 0:
      indent += width
      s = obj.name.decode('utf-8')[:maxkeylen]
      print(s + ' '*(maxkeylen-len(s)) + sep, end='')
    else:
      print(obj.name.decode('utf-8'))

  if len(obj.keys()) > 0:
    for i, key in enumerate(obj.keys()):
      if i>0 or startline:
        print(' '*indent, end='')
      show_uproot_tree(obj[key], indent=indent)
    indent -= width

def processdata(rootfilenames, edfilename, nndatafilename: str, mode: str, 
  nbins = 10, test_ratio: float = 0.2, noisefreqpersqmm: float = 0., noiseTimeRange: float = 5, timeGate = None,
  updateEDF: bool = False, verbose: bool = False, usecupy: bool = False, recorigtime: bool = False, bin3d: bool = False):
  '''
  Чтение и слияние данных из root-файлов в Pandas DataFrame. 
  Транформация данных о событиях в бинированные данные для нейросети и запись их в HDF-файл.
  '''

  # проверка переданных аргументов
  assert(not bin3d and isinstance(nbins, int) and (nbins > 0 and nbins < 100) or 
         bin3d and isinstance(nbins, tuple) and len(nbins)==3 and all([isinstance(n, int) for n in nbins]))
  assert(test_ratio >= 0 and test_ratio <= 1.0)
  assert(noisefreqpersqmm >= 0.)
  assert(noiseTimeRange > 0)
  assert(timeGate is None or timeGate > 0 and timeGate <= noiseTimeRange)
  
  # оценка числа срабатывания для вычитывания HDF
  hitchunksize: int = 10000000

  # диапазон отбора срабатываний по времени
  if timeGate is None:
    timeGate = noiseTimeRange

  # начало отсчета в нс для отбора срабатываний и 3d-бинирования
  timeStart: float = 0
  if recorigtime and timeGate == noiseTimeRange:
    timeStart = -2 # приблизительный сдвиг времени начала срабатываний после восстановления время вылета фотонов

  # полное число бинов
  nTotalBins: int = np.product(nbins) if bin3d else nbins

  # инициализация генератора псевдослучайных числе
  rng = np.random.default_rng(12345)

  def addNoise(partdf: pd.DataFrame, hitdf: pd.DataFrame):
    ''' 
    Генерация темновых срабатываний темнового шума
    '''
    assert(np.isclose(pixel_size*nxpixels_arr+pixel_gap*(nxpixels_arr-1), array_size))
    nevents = partdf.shape[0]

    print(f'    Generate noise with DCR per mm^2 {noisefreqpersqmm}, mean number of hits per event: {munoise:.2f}.')
    startNoise = time.time()

    noisehits = rng.poisson(munoise, partdf.shape[0])
    Ndc = int(noisehits.sum())
    ich = rng.choice(xgrid.size, Ndc)
    xh = xgrid[ich]
    yh = ygrid[ich]
    th = rng.uniform(0, noiseTimeRange, size=Ndc)
    ievent = np.repeat(partdf.index, noisehits)
    ihit = np.zeros(Ndc, dtype=np.int64)
    signalhits = partdf['nhits'].to_numpy()
    zh = hitdf['z_c'].to_numpy()[0]
    index = 0
    for i in range(nevents):
      ihit[index:index+noisehits[i]] = signalhits[i] + np.arange(noisehits[i])
      index += noisehits[i]
    
    noisedf = pd.DataFrame({'x_c': xh, 'y_c': yh, 'z_c': zh, 'time_c': th}, index=pd.MultiIndex.from_arrays((ievent, ihit), names=('entry', 'subentry')))
    hitdf2 = pd.concat((hitdf, noisedf), copy=False).sort_index(level=('entry', 'subentry'))

    print(f'    Noise generation took {time.time()-startNoise:.2f} seconds')

    return hitdf2

  def getAngles(edf: pd.DataFrame):
    '''
    Из координат срабатываний и частиц вычислить углы theta_c и phi_c фотонов и добавить их к edf
    '''
    # начало отсчета для измерения времени преобразования координат в углы
    startAngles = time.time()

    r0 = edf.loc[:, ('x_p', 'y_p', 'z_p')].to_numpy()
    n0 = edf.loc[:, ('nx_p', 'ny_p', 'nz_p')].to_numpy()
    r  = edf.loc[:, ('x_c', 'y_c', 'z_c')].to_numpy()

    # число срабатываний
    N = edf.shape[0]

    if usecupy:
      # Преобразование в CuPy-массивы
      r0, n0, r = cp.asarray(r0), cp.asarray(n0), cp.asarray(r)
      
      # Косинусы и синусы сферических углов направления частицы
      costheta, sintheta = n0[:,2], cp.sqrt(n0[:,0]**2+n0[:,1]**2)
      phi = cp.arctan2(n0[:,1], n0[:,0])
      cosphi, sinphi = cp.cos(phi), cp.sin(phi)

      # Номинальная точка вылета фотонов
      ro = r0 + (W/2+rad_pos)/n0[:,2].reshape(N,1)*n0

      """
      Преобразование в СК частицы
      𝑢𝑥 = cos 𝜃(𝑣𝑥 cos 𝜙 + 𝑣𝑦 sin 𝜙) − 𝑣𝑧 sin 𝜃,
      𝑢𝑦 = −𝑣𝑥 sin 𝜙 + 𝑣𝑦 cos 𝜙,
      𝑢𝑧 = sin 𝜃(𝑣𝑥 cos 𝜙 + 𝑣𝑦 sin 𝜙) + 𝑣𝑧 cos 𝜃.
      """
      # Вектор направления фотона в лабораторной СК
      s = (r-ro)
      snorm = cp.linalg.norm(s, axis=1, keepdims=True)
      v = s / snorm

      # Освобождение памяти
      del r0, n0, ro, r, s

      # Скорректированное время фотона
      if recorigtime:
        edf['time_c'] -= (snorm / speedOfLight_mmperns).reshape(N).get()

      # Матрица преобразования в CK частицы
      U = cp.stack((cp.stack((costheta*cosphi, costheta*sinphi, -sintheta)),
                    cp.stack((-sinphi,         cosphi,          cp.full(N, 0.))),
                    cp.stack((sintheta*cosphi, sintheta*sinphi, costheta)))).transpose(2,0,1)

      # Единичный вектор направления фотона в СК частицы
      u = (U @ v.reshape(N,3,1)).reshape(N,3)

      edf['theta_c'] = cp.arccos(u[:,2]).get()
      edf['phi_c'] = cp.arctan2(-u[:,1], -u[:,0]).get()

    else: # использование NumPy
      # Косинусы и синусы сферических углов направления частицы
      costheta, sintheta = n0[:,2], np.sqrt(n0[:,0]**2+n0[:,1]**2)
      phi = np.arctan2(n0[:,1], n0[:,0])
      cosphi, sinphi = np.cos(phi), np.sin(phi)

      # Номинальная точка вылета фотонов
      ro = r0 + (W/2+rad_pos)/n0[:,2].reshape(N,1)*n0

      """
      Преобразование в СК частицы
      𝑢𝑥 = cos 𝜃(𝑣𝑥 cos 𝜙 + 𝑣𝑦 sin 𝜙) − 𝑣𝑧 sin 𝜃,
      𝑢𝑦 = −𝑣𝑥 sin 𝜙 + 𝑣𝑦 cos 𝜙,
      𝑢𝑧 = sin 𝜃(𝑣𝑥 cos 𝜙 + 𝑣𝑦 sin 𝜙) + 𝑣𝑧 cos 𝜃.
      """

      # Вектор направления фотона в лабораторной СК
      s = (r-ro)
      snorm = np.linalg.norm(s, axis=1, keepdims=True)
      v = s / snorm

      # Освобождение памяти
      #del r0, n0, ro, r, s

      # Скорректированное время фотона
      edf['time_c'] -= (snorm / speedOfLight_mmperns).reshape(N)

      U = np.stack((np.stack((costheta*cosphi, costheta*sinphi, -sintheta)),
                    np.stack((-sinphi,         cosphi,          np.full(N, 0.))),
                    np.stack((sintheta*cosphi, sintheta*sinphi, costheta)))).transpose(2,0,1)

      # Единичный вектор направления фотона в СК частицы
      u = (U @ v.reshape(N,3,1)).reshape(N,3)

      # Сферические углы фотона в СК частицы
      edf['theta_c'] = np.arccos(u[:,2])
      edf['phi_c'] = np.arctan2(-u[:,1], -u[:,0])

      print(f'    Angles evaluated in {time.time()-startAngles:.2f} seconds')

  def formNNdata1d(edf: pd.DataFrame, nndstore: pd.HDFStore, append: bool):
    '''
    Бинирование срабатываний по углу phi_c
    '''

    print(f'  Binning chunk {chunk}')

    # начало отсчета для измерения времени бинирования
    startBinning = time.time()

    # число событий с ненулевым числом хитов
    entries = edf.index.levels[0]
    nentries = entries.size

    # размер бина по phi_c
    assert(isinstance(nbins, int))
    h = 2*np.pi/nbins # шаг

    # отбор срабатываний по углу theta_c
    event_cut = (edf['theta_c'] < np.arccos(1/n_max)*1.2)

    # отбрасывание первого и последнего события в chunk, чтобы избежать частично прочитанных событий в edf
    event_cut = event_cut & (edf.index.get_level_values('entry')!=entries[0]) & (edf.index.get_level_values('entry')!=entries[-1])

    # фильтр срабатываний по времени
    event_cut = event_cut & (edf['time_c'] >= timeStart) & (edf['time_c'] < timeStart+timeGate)

    # выделение информации о частице в отдельный DataFrame
    primdf = edf.loc[event_cut, ('theta_p', 'beta', 'pid', 'momentum')].groupby('entry').mean()

    # вычисление дискретных значений phi_c
    edf['phi_c_bin'] = np.floor((edf['phi_c']+pi)/h)
    edf = edf.loc[event_cut, ('phi_c_bin', 'theta_c', 'time_c')]
    edf.set_index('phi_c_bin', append=True, inplace=True)
    
    # группировка по событиям и дискретным phi_c
    phigroups = edf.groupby(['entry', 'phi_c_bin'])

    # обновление числа событий после наложения условий отбора
    entries = edf.index.levels[0]
    nentries = entries.size
    if verbose:
      print(f'    After cut {nentries} entries left')

    # число срабатываний в бинах
    hitcount = phigroups.count()['theta_c']
    hitcount.name = 'count'

    # среднее в бинах
    meandf = phigroups.mean()
    meandf.rename(columns={'theta_c': 'theta_c_mean', 'time_c': 'time_c_mean'}, inplace=True, errors='raise')
    
    # стандартное отклонение в бинах
    stddf = phigroups.std(0)
    stddf.rename(columns={'theta_c': 'theta_c_std', 'time_c':'time_c_std'}, inplace=True, errors='raise')

    # cоздание dataframe с данными, бинированными по phi_c
    bindf = pd.concat((hitcount,
                       meandf[['theta_c_mean', 'time_c_mean']],
                       stddf[['theta_c_std','time_c_std']]),
                      axis=1)

    # переиндексирование DataFrame с созданием новых строк для бинов без хитов
    key_arrays = [entries.repeat(nbins), np.tile(np.arange(nbins), nentries)]

    # освобождение памяти
    del event_cut, edf, hitcount, meandf, stddf, phigroups

    bindf = bindf.reindex(pd.MultiIndex.from_arrays(key_arrays, names=bindf.index.names), fill_value=0.)
    if verbose:
      print('    Binned data')
      print(bindf)

    if np.isnan(bindf).any().any():
      raise ValueError('There is a NaN value in bindf!')

    # Подготовка данных для нейросети
    data = pd.concat((bindf['count']/50, bindf['theta_c_mean']/pi*2, bindf['theta_c_std']/pi*2, bindf['time_c_mean']/10., bindf['time_c_std']/10.), axis=1)
    
    # Преобразовать индексы по бинам phi_c в колонки
    data = data.unstack(level='phi_c_bin')
    
    if np.isnan(data).any().any():
      raise ValueError('There is a NaN value in NN data!')

    # Проверка совпадения длины DataFrame для частицы и хитов
    assert(primdf.shape[0] == data.shape[0])
    
    data['theta_p_norm'] = primdf['theta_p']/pi*2

    if mode == 'beta':
      labels = pd.Series((primdf['beta']-1/n_max) / (1-1/n_max), index=data.index, name='label', dtype=np.float32)
    elif mode == 'class':
      pids = primdf['pid']
      if upids is None:
        upids = np.unique(pids)
        assert(nhyp == len(upids))
      labels = pd.Series(upids.searchsorted(pids), index=data.index, name='label', dtype=np.int32)
      data['momentum_norm'] = primdf['momentum']/1e4

    if np.isnan(labels).any():
      raise ValueError('There is a NaN value in labels!')

    # преобразование всех данных во float с единичной точностью
    data = data.astype(np.float32)

    print(f'    Binning took {time.time()-startBinning:.2f} seconds')

    if verbose:
      print('    NN data')
      print(data)
      print(labels)
      print(f'    NN data size: {(sys.getsizeof(data)+sys.getsizeof(labels))/1024**2} MiB')
    
    train_data = test_data = train_labels = test_labels = None
    if test_ratio>0.0 and test_ratio<1.0:
      print(f'    Splitting data into train and test sets with test size ratio {test_ratio}')
      train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=test_ratio)
      if verbose:
        print(f'    Shapes of train and test data and labels:', train_data.shape, train_labels.shape, test_data.shape, test_labels.shape)
    elif test_ratio==0.0:
      train_data = data
      train_labels = labels
    elif test_ratio==1.0:
      test_data = data
      test_labels = labels

    print(f'    Saving NN data chunk {chunk}...')
    #hdfmode = ('a' if append else 'w')
    # Сохранение данных для нейросети
    if train_data is not None:
      nndstore.put('train_data', train_data, format='table', append=append)
      nndstore.put('train_labels', train_labels, format='table', append=append)
    if test_data is not None:
      nndstore.put('test_data', test_data, format='table', append=append)
      nndstore.put('test_labels', test_labels, format='table', append=append)

  def formNNdata3d(edf: pd.DataFrame, nndstore: pd.HDFStore, append: bool):
    '''
    Бинирование срабатываний по углам phi_c, theta_c и времени вылета фотона
    '''

    print(f'  Binning chunk {chunk}')

    # число бинов
    assert(isinstance(nbins, tuple) and len(nbins)==3)
    nPhiBins, nCosThetaBins, nTimeBins = nbins

    # начало отсчета для измерения времени бинирования
    startBinning = time.time() #0

    # число событий с ненулевым числом хитов
    entries = edf.index.levels[0]
    nentries = entries.size

    # размер бина по phi_c
    hphi = 2*np.pi/nPhiBins # шаг
    
    # размер бина по theta_c
    maxTheta = np.arccos(1/n_max)*1.2
    hctheta = (1-np.cos(maxTheta)) / nCosThetaBins

    # размер бина по time_c
    htime = timeGate / nTimeBins

    # отбор срабатываний по углу theta_c
    event_cut = (edf['theta_c'] < maxTheta)

    # отбрасывание первого и последнего события в chunk, чтобы избежать частично прочитанных событий в edf
    event_cut = event_cut & (edf.index.get_level_values('entry')!=entries[0]) & (edf.index.get_level_values('entry')!=entries[-1])

    # фильтр срабатываний по времени
    event_cut = event_cut & (edf['time_c'] >= timeStart) & (edf['time_c'] < timeStart+timeGate)

    edf = edf[event_cut]

    eventCutTime = time.time() #1

    # выделение информации о частице в отдельный DataFrame
    primdf = edf.loc[:, ('theta_p', 'beta', 'pid', 'momentum')].groupby('entry').mean()

    # вычисление дискретных значений бинов
    edf['bin'] = nCosThetaBins * nTimeBins * np.minimum(np.floor((edf.loc[:,'phi_c']+pi) / hphi), nPhiBins-1) + \
                 nTimeBins * np.minimum(np.floor((1-np.cos(edf.loc[:, 'theta_c'])) / hctheta), nCosThetaBins-1) + \
                 np.minimum(np.maximum(np.floor((edf.loc[:, 'time_c']-timeStart) / htime), 0), nTimeBins)
    edf.reset_index('subentry', inplace=True)
    edf = edf.loc[:, ('bin', 'subentry')]
    edf.set_index('bin', append=True, inplace=True)

    binEdgeTime = time.time() #2

    # обновление числа событий после наложения условий отбора
    entries = edf.index.levels[0]
    nentries = entries.size
    if verbose:
      print(f'    After cut {nentries} entries left')

    # число срабатываний в бинах
    bindf = edf.groupby(['entry', 'bin']).count().astype(np.float32) / 10
    bindf.columns = ['count']

    groupCountsTime = time.time() #3

    # освободить память
    del event_cut, edf

    # переиндексирование DataFrame с созданием новых строк для бинов без хитов
    key_arrays = [entries.repeat(nTotalBins), np.tile(np.arange(nTotalBins), nentries)]

    bindf = bindf.reindex(pd.MultiIndex.from_arrays(key_arrays, names=bindf.index.names), fill_value=0.)

    reindexTime = time.time() #4

    if np.isnan(bindf).any().any():
      raise ValueError('There is a NaN value in bindf!')

    # Подготовка данных для нейросети
    data = bindf
    
    # Преобразовать индексы по бинам в колонки
    data = data.unstack(level='bin')

    unstackTime = time.time() #5
    
    if np.isnan(data).any().any():
      raise ValueError('There is a NaN value in NN data!')

    # Проверка совпадения длины DataFrame для частицы и хитов
    assert(primdf.shape[0] == data.shape[0])
    
    data['theta_p_norm'] = (primdf['theta_p']/pi*2).astype(np.float32)

    if mode == 'beta':
      labels = pd.Series((primdf['beta']-1/n_max) / (1-1/n_max), index=data.index, name='label', dtype=np.float32)
    elif mode == 'class':
      pids = primdf['pid']
      if upids is None:
        upids = np.unique(pids)
        assert(nhyp == len(upids))
      labels = pd.Series(upids.searchsorted(pids), index=data.index, name='label', dtype=np.int32)
      data['momentum_norm'] = (primdf['momentum']/1e4).astype(np.float32)

    if np.isnan(labels).any():
      raise ValueError('There is a NaN value in labels!')

    endBinning = time.time() #6

    print(f'    Binning took {endBinning-startBinning:.2f} seconds of which')
    print(f'      Event cut: {eventCutTime-startBinning:.2f} seconds')
    print(f'      Bin edges: {binEdgeTime-eventCutTime:.2f} seconds')
    print(f'      Hit counting: {groupCountsTime-binEdgeTime:.2f} seconds')
    print(f'      Reindexing: {reindexTime-groupCountsTime:.2f} seconds')
    print(f'      Unstacking: {unstackTime-reindexTime:.2f} seconds')
    print(f'      Adding particle data: {endBinning-unstackTime:.2f} seconds')

    if verbose:
      print('    NN data')
      print(data)
      print(labels)
      print(f'    NN data size: {(sys.getsizeof(data)+sys.getsizeof(labels))/1024**2} MiB')
    
    train_data = test_data = train_labels = test_labels = None
    if test_ratio>0.0 and test_ratio<1.0:
      print(f'    Splitting data into train and test sets with test size ratio {test_ratio}')
      train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=test_ratio)
      if verbose:
        print(f'    Shapes of train and test data and labels:', train_data.shape, train_labels.shape, test_data.shape, test_labels.shape)
    elif test_ratio==0.0:
      train_data = data
      train_labels = labels
    elif test_ratio==1.0:
      test_data = data
      test_labels = labels

    print(f'    Saving NN data chunk {chunk}...')
    # Сохранение данных для нейросети
    if train_data is not None:
      nndstore.put('train_data', train_data, format='table', append=append)
      nndstore.put('train_labels', train_labels, format='table', append=append)
    if test_data is not None:
      nndstore.put('test_data', test_data, format='table', append=append)
      nndstore.put('test_labels', test_labels, format='table', append=append)

  def readInfoFromRoot():
    '''
    Получение информации о моделировании из ROOT-файла в виде idf data frame
    '''
    # Объявление названий колонок данных и определение колонок для записи 
    idf_rename_map = {'m_num_events': 'nevents',
                      'm_z_dis': 'zdis',
                      'm_layers': 'nlayers',
                      'm_size': 'array_size', 'm_gap': 'array_gap',
                      'm_chip_size': 'pixel_size', 'm_chip_pitch': 'pixel_gap',
                      'm_chip_num_size': 'pixel_numx',
                      'm_num_side_x': 'nxarrays', 'm_num_side_y': 'nyarrays',
                      'm_focal_length': 'distance',
                      'm_trg_window': 'trg_window_ns',
                      'W': 'W', 'n_mean': 'n_mean', 'n_max': 'n_max', 'nhyp': 'nhyp'}

    # Открытие ROOT-файла с данными используя Uproot https://github.com/scikit-hep/uproot3
    with uproot.open(rootfilenames[0]) as file:
      idf = file[b'info_sim'].pandas.df('*', flatten=False)

    idf.rename(columns=idf_rename_map, inplace=True, errors='ignore')

    # Получение параметров (многослойного) радиатора одинаковых для всех файлов
    n_l = idf.at[0,'m_layers.first'] #показатели преломления слоёв
    w_l = idf.at[0,'m_layers.second'] #толщина слоёв радиатора

    W = w_l.sum() # суммарная толщина всех слоёв
    n_mean = n_l.mean() # средний показатель преломления
    n_max = n_l.max() # максимальный показатель преломления

    idf['W'] = W
    idf['n_mean'] = n_mean
    idf['n_max'] = n_max
    idf['nhyp'] = 1

    idf['nevents'] = 0

    idf = idf[idf_rename_map.values()]

    return idf

  def readInfoFromHDF():
    idf = pd.read_hdf(edfilename, 'idf')
    return idf

  def genChunkFromRoot():
    '''
    Генерация событий из ROOT-файла
    '''
    part_rename_map = {'m_hits': 'nhits', 'm_pos_primary._0': 'x_p', 'm_pos_primary._1': 'y_p', 'm_pos_primary._2': 'z_p', 
                       'm_dir_primary._0': 'nx_p', 'm_dir_primary._1': 'ny_p', 'm_dir_primary._2': 'nz_p', 
                       'm_id_primary': 'pid','m_beta_primary': 'beta', 'm_theta_primary': 'theta_p', 'm_phi_primary': 'phi_p',
                       'm_momentum_primary': 'momentum'}
    hit_rename_map = {'m_hits.m_photon_pos_chip._0': 'x_c', 'm_hits.m_photon_pos_chip._1': 'y_c', 'm_hits.m_photon_pos_chip._2': 'z_c',
                      'm_hits.m_photon_time': 'time_c'}
    #edfcolstosave = ['theta_p', 'phi_p', 'beta', 'pid', 'momentum', 'mass', 'x_c', 'y_c', 'z_c', 'theta_c', 'phi_c', 'time_c'] # for debugging
    edfcolstosave = ['theta_p', 'phi_p', 'beta', 'pid', 'momentum', 'mass', 'theta_c', 'phi_c', 'time_c']

    print(f'Reading and processing data from ROOT-files...')
    print(f'Event chunk size {eventchunksize}')

    entryShift = 0

    # Цикл по ROOT-файлам
    for fn in rootfilenames:
      with uproot.open(fn) as file:
        nFileEvents = file[b'info_sim/info_gen/m_num_events'].array()[0]

      print(f'Processing ROOT file {fn} with {nFileEvents} simulated events...', flush=True)
      # Цикл по кускам в ROOT-файле
      for partdf, hitdf in zip(uproot.pandas.iterate(fn, "raw_data", part_rename_map.keys(), entrysteps=eventchunksize), 
                               uproot.pandas.iterate(fn, "raw_data", hit_rename_map.keys(), entrysteps=eventchunksize, flatten=True)):
        print(f'\n  Processing chunk {chunk}...')

        # Переименование колонок
        partdf.rename(columns=part_rename_map, inplace=True, errors='raise')
        hitdf.rename(columns=hit_rename_map, inplace=True, errors='raise')

        # Генерация и добавление шумовых срабатываний
        if noisefreqpersqmm > 0:
          hitdf = addNoise(partdf, hitdf)

        print(f'    {hitdf.index.levels[0].size} entries with {hitdf.shape[0]} hits to process')

        # Вычисление массы частиц
        pdg = partdf['pid'].to_numpy()
        used_pdgids = np.unique(pdg)
        nhyp = used_pdgids.size # число гипотез (разных типов частиц)
        idf['nhyp'] = max(nhyp, idf.at[0, 'nhyp'])

        particles = dict([(pdgid, Particle.from_pdgid(pdgid)) for pdgid in used_pdgids])
        partdf['mass'] = [particles[code].mass for code in pdg]

        if verbose:
          for pdgid, particle in particles.items():
            print(f'    PDGCODE={pdgid}: {particle.name}')

        # Слияние данных событий и срабатываний
        edf = hitdf.join(partdf, on='entry')

        # Число событий с ненулевым числом хитов
        idf['nevents'] += edf.index.levels[0].size

        # Вычисление углов theta_c и phi_c, коррекция времени
        getAngles(edf)

        edf['time_c']

        # Переиндексация для создания уникальных номеров событий
        if entryShift > 0:
          # Номера событий и хитов текущего DataFrame
          entries = edf.index.get_level_values('entry') + entryShift
          subentries = edf.index.get_level_values('subentry')
          edf.set_index(pd.MultiIndex.from_arrays((entries, subentries), names=edf.index.names), inplace=True)

        edf = edf[edfcolstosave]

        if verbose:
          print(edf)

        if edfstore is not None:
          print(f'    Saving EDF chunk...')
          edfstore.put('edf', edf, format='table', append=(chunk>0))

        yield edf

      # сдвиг начального номера события для следующего файла
      entryShift += nFileEvents
      
  def genChunkFromHDF():
    edfstore = pd.HDFStore(edfilename, 'r')
    print(f'Reading data from HDF {edfilename}...')
    print('HDF content:')
    print(edfstore.info())
    print(f'Hit chunk size {hitchunksize}')

    for edf in edfstore.select('edf', chunksize=hitchunksize):
      print(f'\n  Processing chunk {chunk}...')
      print(f'    {edf.index.levels[0].size} entries with {edf.shape[0]} hits to process')
      print(f'    EDF size: {edf.size}')
      yield edf
    
  #==========================================================================================
  print('='*100)

  genEDF = None
  edfExists = edfilename is not None and os.access(edfilename, os.R_OK)
  updateEDF = edfilename is not None and updateEDF
  edfstore = None
  if not edfExists or updateEDF:
    print(f'Read simulation data from {rootfilenames}.')
    if isinstance(rootfilenames, str):
      rootfilenames = [rootfilenames]

    if updateEDF:
      if os.access(edfilename, os.R_OK):
        os.remove(edfilename)

      print(f'Write edf data frame to {edfilename}.')
      edfstore = pd.HDFStore(edfilename, 'w')

    genEDF = genChunkFromRoot

    idf = readInfoFromRoot()

  elif edfExists:
    print(f'Read EDF from {edfilename}')

    genEDF = genChunkFromHDF

    idf = readInfoFromHDF()

  # определения параметра фотодетектора для генерации темнового шума
  pixel_size, pixel_gap = idf.at[0,'pixel_size'], idf.at[0, 'pixel_gap']
  array_size, array_gap = idf.at[0,'array_size'], idf.at[0,'array_gap']
  nxpixels_arr = idf.at[0, 'pixel_numx']
  nxpixels_tot = idf.at[0,'nxarrays']*nxpixels_arr
  igrid = np.arange(nxpixels_tot//2)
  xpnts = array_gap/2 + (igrid//nxpixels_arr)*(array_size+array_gap) + (igrid%nxpixels_arr)*(pixel_size+pixel_gap) + pixel_size/2
  xpnts = np.sort(np.append(-xpnts, xpnts))
  xgrid, ygrid = np.meshgrid(xpnts, xpnts)
  xgrid = xgrid.reshape(xgrid.size)
  ygrid = ygrid.reshape(ygrid.size)

  # максимальный показатель преломления
  n_max = float(idf['n_max'])

  # толщина радиатора
  W = float(idf['W'])

  # положение радиатора
  rad_pos = float(idf['zdis'])

  # оценка числа событий для вычитывания из ROOT-файла
  munoise = noiseTimeRange*1e-9*noisefreqpersqmm*(pixel_size**2)*(nxpixels_tot**2)
  eventchunksize: int = int(min(hitchunksize/(50+munoise), 25000000/nTotalBins))
  
  print(f'Write NN data to {nndatafilename}.')
  nndstore = pd.HDFStore(nndatafilename, 'w')

  if usecupy:
    print('Use CuPy for computations on GPU')

  # выбор метода бинирования
  if bin3d:
    print(f'Bin data on phi_c ({nbins[0]} bins), theta_c ({nbins[1]} bins) and time_c ({nbins[2]} bins)')
    formNNdata = formNNdata3d
  else:
    print(f'Bin data on phi_c only ({nbins} bins)')
    formNNdata = formNNdata1d

  startJob = time.time()

  try:
    chunk = 0

    for edf in genEDF():
      nhyp = float(idf['nhyp'])

      formNNdata(edf, nndstore, append=(chunk>0))

      chunk += 1

  finally:
    if edfstore is not None:
      edfstore.close()
    nndstore.close()

    if updateEDF:
      # сохранение idf вместе с edf
      print(f'Storing idf to {edfilename}')
      idf.to_hdf(edfilename, 'idf', mode='a', append=False)

    # сохранение idf вместе с NN data
    print(f'Storing idf to {nndatafilename}')
    idf.to_hdf(nndatafilename, 'idf', mode='a', append=False)

  print(f'Job finished in {time.time()-startJob:.2f} seconds')

  print('='*100)

def plothits(fn, events=10000, hits=500000):
  '''
  Отрисовка срабатываний в ROOT-файле или DataFrame
  '''  
  figxy = plt.figure(figsize=(10,10))
  if fn[-5:]=='.root':
    with uproot.open(fn) as f:
      xh = f[b'raw_data/event/m_hits/m_hits.m_photon_pos_chip._0'].array(entrystop=events).flatten()
      yh = f[b'raw_data/event/m_hits/m_hits.m_photon_pos_chip._1'].array(entrystop=events).flatten()
      zh = f[b'raw_data/event/m_hits/m_hits.m_photon_pos_chip._2'].array(entrystop=events).flatten()
      th = f[b'raw_data/event/m_hits/m_hits.m_photon_time'].array(entrystop=events).flatten()
      chip_id = f[b'raw_data/event/m_hits/m_hits.m_photon_id_chip'].array(entrystop=events).flatten()
      pmt_id = f[b'raw_data/event/m_hits/m_hits.m_photon_id_pmt'].array(entrystop=events).flatten()
      for i in range(20):
        print('xh=%7.2f yh=%7.2f zh=%7.2f th=%7.2f chipID=%2d pmtID=%3d' % (xh[i], yh[i], zh[i], th[i], chip_id[i], pmt_id[i]))
      npmts = max(pmt_id)+1
      cmap = plt.get_cmap('viridis')
      for a in range(npmts):
        mask = (pmt_id==a)
        if not mask.any():
          continue
        plt.plot(xh[mask], yh[mask],'r.', ms=1, color=cmap(float(a)/npmts))
  elif fn[-3:]=='.h5':
    edf = next(iter(pd.read_hdf(fn, 'edf', 'r', chunksize=hits)))
    th = edf['time_c'].to_numpy()
    theta_c = edf['theta_c'].to_numpy()
    if 'x_c' in edf.columns:
      xh = edf['x_c'].to_numpy()
      yh = edf['y_c'].to_numpy()
      zh = edf['z_c'].to_numpy()
      for i in range(20):
        print('xh=%7.2f yh=%7.2f zh=%7.2f th=%7.2f' % (xh[i], yh[i], zh[i], th[i]))
      plt.plot(xh, yh,'b.', ms=0.5)
  else:
    print(f'Unknown file extension: {fn}')
    return

  if 'xh' in locals():
    print("(Xmin, Xmax):", xh.min(), xh.max())
    print("(Tmin, Tmax):", th.min(), th.max())
    plt.plot([-400, 400], [0, 0], 'k-')
    plt.plot([0, 0], [-400, 400], 'k-')
    plt.xlim(-30, 30)
    plt.ylim(-30, 30)
    plt.xlabel('Xhit, mm')
    plt.ylabel('Yhit, mm')

  figtime = plt.figure(figsize=(10,6))
  plt.hist(th, bins=100)
  plt.xlabel('Thit, ns', fontsize=16)
  plt.ylabel('Hits/0.05ns', fontsize=16)
  plt.semilogy()

  figtheta = None
  if 'theta_c' in locals():
    figtheta = plt.figure(figsize=(10,6))
    plt.hist(edf['theta_c'], bins=100)
    plt.xlabel('$theta_c$', fontsize=16)
    plt.ylabel('Hits', fontsize=16)
    plt.semilogy()
  
  return figxy, figtime, figtheta

def get_sim_pars(filename: str):
  '''
  Чтение параметров моделирования n_max и nhyp из idf
  '''
  idf = pd.read_hdf(filename, 'idf')
  return (float(idf['n_max']), int(idf['nhyp']))

nndataRepo = dict()

def get_nndata_format(nndatafilename: str):
  '''
  Получение типа и размерности входных данных нейросети из сохраненного DataFrame
  '''
  global nndataRepo
  if nndatafilename in nndataRepo:
    _, data, _, labels = nndataRepo[nndatafilename]
  else:
    try:
      data = pd.read_hdf(nndatafilename, 'test_data', stop=1)
      labels = pd.read_hdf(nndatafilename, 'test_labels', stop=1)
    except KeyError:
      data = pd.read_hdf(nndatafilename, 'train_data', stop=1)
      labels = pd.read_hdf(nndatafilename, 'train_labels', stop=1)

  dtypes = (data.dtypes[0], labels.dtype)
  shapes = (data.shape[1:], ())

  return dtypes, shapes

def load_nndata(nndatafilename: str, types: str = 'all', entries: int = None):
  '''
  Загрузка входных данных нейросети из HDF-файла. Использование ранее загруженных данных при наличии.
    nndatafilename: имя HDF-файла
    types: одно из 'all', 'train', 'test'
  '''
  global nndataRepo

  assert(types in ('all', 'train', 'test'))

  test_data = test_labels = train_data = train_labels = None
  
  if nndatafilename in nndataRepo:
    print(f'Using previously loaded NN data from {nndatafilename}.')
    train_data, test_data, train_labels, test_labels = nndataRepo[nndatafilename]

  if (test_data is None or test_labels is None or 
      entries is not None and test_data.shape[0]<entries) and (types=='all' or types=='test'):
    print(f'Loading test NN data from {nndatafilename}...')
    try:
      test_data = pd.read_hdf(nndatafilename, 'test_data', stop=entries)
      test_labels = pd.read_hdf(nndatafilename, 'test_labels', stop=entries)
    except KeyError:
      print('No test_data or test_labels found!')
    else:
      nndataRepo[nndatafilename] = (train_data, test_data, train_labels, test_labels)

  if (train_data is None or train_labels is None or
      entries is not None and train_data.shape[0]<entries) and (types=='all' or types=='train'):
    print(f'Loading train NN data from {nndatafilename}...')
    try:
      train_data = pd.read_hdf(nndatafilename, 'train_data', stop=entries)
      train_labels = pd.read_hdf(nndatafilename, 'train_labels', stop=entries)
    except KeyError:
      print('No train_data or train_labels found!')
    else:
      nndataRepo[nndatafilename] = (train_data, test_data, train_labels, test_labels)

  if types=='all':
    return train_data, test_data, train_labels, test_labels
  elif types=='train':
    return train_data, train_labels
  else:
    return test_data, test_labels

def build_model(input_shape, mode, hlayers=(10,10), activation='relu', nhyp=2, dropout=0.):
    '''
    Построение модели нейросети
    '''
    assert(len(hlayers)>0 and all([type(l) is int for l in hlayers]))
    model = tf.keras.Sequential()
    print(f'Build TF model with hidden layer configuration {hlayers} and input shape {input_shape}')
    
    if dropout > 0. and dropout < 1.:
      model.add(tf.keras.layers.Dropout(dropout, noise_shape=input_shape, input_shape=input_shape))
      model.add(tf.keras.layers.Dense(hlayers[0], activation=activation))
    else:
      model.add(tf.keras.layers.Dense(hlayers[0], activation=activation, input_shape=input_shape))
    for hl in hlayers[1:]:
      model.add(tf.keras.layers.Dense(hl, activation=activation))
    if mode=='class':
      model.add(tf.keras.layers.Dense(nhyp, activation='softmax'))
    else: # mode == 'beta'
      model.add(tf.keras.layers.Dense(1))
    
    optimizer = tf.keras.optimizers.Adam()
    
    if mode == 'class':
      model.compile(optimizer=optimizer, #это то, как модель обновляется на основе данных, которые она видит, и функции потери
                    loss='sparse_categorical_crossentropy', #измеряет насколько точная модель во время обучения
                    metrics=['accuracy']) #используется для контроля за этапами обучения и тестирования
    elif mode == 'beta':
      model.compile(optimizer=optimizer, #это то, как модель обновляется на основе данных, которые она видит, и функции потери
                    loss='mse', #измеряет насколько точная модель во время обучения
                    metrics=['mse']) #используется для контроля за этапами обучения и тестирования

    return model

def fitGaus(data):
  '''
  Подгонка распределения остатков скорости гауссом и возврат сигмы.
  '''
  #alpha = 0.01
  #q1, q2 = list(data.quantile([alpha,1-alpha]))
  cutdata = data
  for i in range(2):
    mean, sd = cutdata.mean(), cutdata.std()
    q1, q2 = mean-3*sd, mean+3*sd
    cutdata = cutdata[(q1<data)&(data<q2)]

  blh = BinnedLH(gaussian, cutdata, bins=100)
  m = Minuit(blh, mean=cutdata.mean(), sigma=cutdata.std())
  m.errordef = Minuit.LIKELIHOOD

  m.migrad() # fit
  sigma = m.values[1]

  return sigma

def gausSigmaWithError(data, iname: str, onames: tuple):
  '''
  Подгонка распределения остатков скорости гауссом и возврат сигмы с ошибкой.
  '''
  assert(len(onames)==3)
  series = data[iname]
  allEntries = series.size
  for i in range(2):
    mean, sd = series.mean(), series.std()
    q1, q2 = mean-3*sd, mean+3*sd
    series = series[(q1<series)&(series<q2)]
  cutoutfrac = 1-series.size/allEntries

  blh = BinnedLH(gaussian, series, bins=100)
  m = Minuit(blh, mean=series.mean(), sigma=series.std())
  m.errordef = Minuit.LIKELIHOOD

  m.migrad() # fit

  return pd.Series((m.values[1], m.errors[1], cutoutfrac), index=onames)

def evaluate_model(test_data, test_labels, model, mode, transfunc, outdir, show=False, savefig=False):
  '''
  Функция оценки работы нейросети на тестовых данных с получением точности восстановленных параметров.
  '''
  prediction = None
  print(f'Evaluate model on {test_data.shape[0]} events')
  t0 = time.perf_counter()
  loss, metric = model.evaluate(test_data, test_labels, verbose=2)
  print(f'Evaluation took {time.perf_counter()-t0:.2f} seconds')

  t0 = time.perf_counter()
  if mode=='beta':
    assert(callable(transfunc))
    
    prediction = model.predict(test_data).flatten()
    elapsed = time.perf_counter()-t0
    print(f'Prediction took {elapsed:.2f} seconds')

    betas = transfunc(test_labels.values)
    pbetas = transfunc(prediction)
    beta_residual = pbetas - betas
    print(f'Systematic shift of beta w.r.t. MC truth: {beta_residual.mean():.3g}')
    print(f'Standard deviation of beta residuals: {beta_residual.std():.3g}')

    if show:
      # Разбиение данных по бинам по скорости и углу падения
      thetas = test_data['theta_p_norm'].values.flatten()*90 # градусы
      beta_min, beta_max = betas.min(), betas.max()+1e-5
      theta_min, theta_max = thetas.min(), thetas.max()+1e-4
      beta_set, theta_set = np.unique(betas), np.unique(thetas)
      if beta_set.size < 20:
        nbetapnts = beta_set.size
        beta_pnt = betas
        dbeta = 0
      else:
        nbetapnts = 10
        beta_pnt = beta_min + (beta_max-beta_min) * (np.floor(nbetapnts * (betas-beta_min) / (beta_max-beta_min)) + 0.5) / nbetapnts
        beta_set = np.unique(beta_pnt)
        dbeta = (beta_max - beta_min) / nbetapnts / 2

      if theta_set.size < 20:
        nthetapnts = theta_set.size
        theta_pnt = thetas
        dtheta = 0
      else:
        nthetapnts = 9
        theta_pnt = theta_min + (theta_max-theta_min) * (np.floor(nthetapnts * (thetas-theta_min) / (theta_max-theta_min)) + 0.5) / nthetapnts
        theta_set = np.unique(theta_pnt)
        dtheta = (theta_max - theta_min) / nthetapnts / 2

      resdf = pd.DataFrame({'beta_pnt': beta_pnt, 'theta_pnt': theta_pnt, 'residual': beta_residual}, dtype=np.float64)
      resdfgrp = resdf.groupby(['beta_pnt', 'theta_pnt'])
      accdf = resdfgrp.mean()
      accdf.columns = ['shift']
      accdf['rms'] = resdfgrp.std()
      print('Minimum standard deviation among theta_p-beta bins:', min(accdf['rms']))
      
      accdf = accdf.join(resdfgrp.apply(gausSigmaWithError, iname='residual', onames=('sigma', 'sigma_error', 'cutout_fraction')))
      print('Minimum fitted gaussian sigma among theta_p-beta bins:', min(accdf['sigma']))

      accdf.reset_index(('beta_pnt', 'theta_pnt'), inplace=True)
      accdf.sort_values(['theta_pnt', 'beta_pnt'], ascending=True, inplace=True)

      # Нарисовать распределения для оценки точности восстановления скорости
      plt.figure(figsize=(10,10))
      plt.subplot(221)
      plt.scatter(betas, pbetas, 0.5, 'b')
      plt.xlabel(r'MC-truth $\beta$', fontsize=16)
      plt.ylabel(r'Predicted $\beta$', fontsize=16)

      plt.subplot(222)
      plt.hist(beta_residual, bins=100, color='b', log=True)
      plt.xlabel(r'$\beta$ residuals', fontsize=16)
      plt.ylabel(r'Events', fontsize=16)

      plt.subplot(223)
      plt.scatter(betas, beta_residual, 0.5, 'b')
      plt.xlabel(r'$\beta$', fontsize=16)
      plt.ylabel(r'$\beta$ residuals', fontsize=16)

      plt.subplot(224)
      plt.scatter(thetas, beta_residual, 0.5, 'r')
      plt.xlabel(r'$\theta_p, deg$', fontsize=16)
      plt.ylabel(r'$\beta$ residuals', fontsize=16)
      plt.tight_layout()
      if savefig:
        plt.savefig(outdir + os.sep + mode + '_residuals.png')

      plt.figure(figsize=(10,10))
      rmsbeta_map = np.reshape(accdf['rms'].to_numpy(), (nthetapnts, nbetapnts))
      plt.imshow(rmsbeta_map,
                 extent=[beta_min, beta_max, theta_min, theta_max],
                 aspect='auto', 
                 origin='lower', # default origin of the plot is the upper-left
                 cmap='viridis')
      def tickpos(vmin, vmax, num):
        hstep = (vmax-vmin)/num/2
        return np.linspace(vmin+hstep, vmax-hstep, num)
      if dbeta == 0.:
        plt.xticks(ticks=tickpos(beta_min, beta_max, nbetapnts), labels=[f'{b:.3f}' for b in beta_set])
      if dtheta == 0.:
        plt.yticks(ticks=tickpos(theta_min, theta_max, nthetapnts), labels=[f'{t:.1f}' for t in theta_set])
      plt.colorbar()
      plt.xlabel(r'$\beta$', fontsize=16)
      plt.ylabel(r'$\theta_p, deg$', fontsize=16)
      plt.title(r'Standard deviation of $\beta$ residuals', fontsize=20)
      plt.tight_layout()
      if savefig:
        plt.savefig(outdir + os.sep + mode+'_sd_resolution_2d.png')

      plt.figure(figsize=(10,10))
      sigbeta_map = np.reshape(accdf['sigma'].to_numpy(), (nthetapnts, nbetapnts))
      plt.imshow(sigbeta_map,
                 extent=[beta_min, beta_max, theta_min, theta_max],
                 aspect='auto', 
                 origin='lower', # default origin of the plot is the upper-left
                 cmap='viridis')
      if dbeta == 0.:
        plt.xticks(ticks=tickpos(beta_min, beta_max, nbetapnts), labels=[f'{b:.3f}' for b in beta_set])
      if dtheta == 0.:
        plt.yticks(ticks=tickpos(theta_min, theta_max, nthetapnts), labels=[f'{t:.1f}' for t in theta_set])
      plt.colorbar()
      plt.xlabel(r'$\beta$', fontsize=16)
      plt.ylabel(r'$\theta_p, deg$', fontsize=16)
      plt.title(r'Sigma of gaussian fitted to truncated $\beta$ residuals', fontsize=20)
      plt.tight_layout()
      if savefig:
        plt.savefig(outdir + os.sep + mode+'_gaus_resolution_2d.png')
      
      plt.figure(figsize=(10,10))
      cutoutfrac_map = np.reshape(accdf['cutout_fraction'].to_numpy(), (nthetapnts, nbetapnts))
      plt.imshow(cutoutfrac_map,
                 extent=[beta_min, beta_max, theta_min, theta_max],
                 aspect='auto', 
                 origin='lower', # default origin of the plot is the upper-left
                 cmap='viridis')
      if dbeta == 0.:
        plt.xticks(ticks=tickpos(beta_min, beta_max, nbetapnts), labels=[f'{b:.3f}' for b in beta_set])
      if dtheta == 0.:
        plt.yticks(ticks=tickpos(theta_min, theta_max, nthetapnts), labels=[f'{t:.1f}' for t in theta_set])
      plt.colorbar()
      plt.xlabel(r'$\beta$', fontsize=16)
      plt.ylabel(r'$\theta_p, deg$', fontsize=16)
      plt.title(r'Cut-out fraction after truncation of $\pm 3\sigma$', fontsize=20)
      plt.tight_layout()
      if savefig:
        plt.savefig(outdir + os.sep + mode+'_cutout_fraction_2d.png')
      
      plt.figure(figsize=(10,8))
      sbmax = 0      
      for theta in theta_set:
        mask = np.isclose(accdf['theta_pnt'], theta)
        b, sb, esb = accdf.loc[mask, 'beta_pnt'].to_numpy(), accdf.loc[mask, 'sigma'].to_numpy(), accdf.loc[mask, 'sigma_error'].to_numpy()
        if sbmax < sb.max():
          sbmax = sb.max()
        plt.errorbar(b, sb, xerr=dbeta, yerr=esb, fmt='-o', label=rf'$\theta_p =$ ({theta-dtheta:.1f}..{theta+dtheta:.1f}$^0$)')
      plt.xlabel(r'$\beta$', fontsize=16)
      plt.ylabel(r'$\sigma_{\beta}$', fontsize=16)
      plt.title(r'Gaussian velocity resolution for different incident angles', fontsize=18)
      plt.legend(fontsize=16)
      plt.ylim(0, 1.2*sbmax)
      plt.grid()
      plt.tight_layout()
      if savefig:
        plt.savefig(outdir + os.sep + mode+'_gaus_resolution_1d.png')

  elif mode=='class':
    prediction = model.predict(test_data)
    nhyp = prediction.shape[1]
    prediction = prediction.argmax(axis=1).flatten()
    elapsed = time.perf_counter()-t0
    print(f'Prediction took {elapsed:.2f} seconds')

    true_pid_rates = np.zeros(nhyp)
    for ihyp in range(nhyp):
      hyp_mask = (test_labels==ihyp)
      true_pid_rates[ihyp] = ((prediction[hyp_mask]==ihyp).sum()/hyp_mask.sum())
    plt.figure(figsize=(10,8))
    plt.bar(range(nhyp), true_pid_rates, tick_label=[r'$\pi$', r'$\mu$'])
    plt.ylim(true_pid_rates.min()-(1-true_pid_rates.min())*0.1,1)
    plt.title('Particle true identification probability')
  
  return metric, prediction