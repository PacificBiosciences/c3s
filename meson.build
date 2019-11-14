project(
  'c3s',
  ['c','cpp'],
  version : '0.0.1',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++14',
    'c_std=gnu11'],
  license : 'BSD-3',
  meson_version : '>= 0.46.0')

# CXXFLAGS
c3s_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-unused-command-line-argument',
  '-Wno-missing-field-initializers',
  '-Wno-unused-const-variable',
  '-Wno-c99-extensions']
    if cpp.has_argument(cflag)
      c3s_flags += cflag
    endif
endforeach

c3s_c_flags = []
c = meson.get_compiler('c')
foreach cflag: [
  '-Wc++-compat',
  '-Wno-unused-parameter',
  '-Wno-unused-variable',
  '-Wno-sign-compare']
    if c.has_argument(cflag)
      c3s_c_flags += cflag
    endif
endforeach

if get_option('sse41')
  if cpp.has_argument('-msse4.1')
    c3s_flags += '-msse4.1'
  else
    error('Requested SSE4.1, but the C++ compiler does not support -msse4.1')
  endif

  if c.has_argument('-msse4.1')
    c3s_c_flags += '-msse4.1'
  else
    error('Requested SSE4.1, but the C compiler does not support -msse4.1')
  endif
endif

# dependencies
## threads
c3s_thread_dep = dependency('threads', required : true)
## boost
c3s_boost_dep = dependency('boost', required : true)
## zlib
c3s_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])
## htslib
c3s_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])
## pbbam
c3s_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])
## pbcopper
c3s_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])

c3s_lib_deps = [
  c3s_thread_dep,
  c3s_boost_dep,
  c3s_zlib_dep,
  c3s_htslib_dep,
  c3s_pbbam_dep,
  c3s_pbcopper_dep
]

# sources
c3s_third_party_cpp_sources = files([
  'src/third-party/spoa/alignment_engine.cpp',
  'src/third-party/spoa/graph.cpp',
  'src/third-party/spoa/sequence.cpp',
  'src/third-party/spoa/simd_alignment_engine.cpp',
  'src/third-party/spoa/sisd_alignment_engine.cpp',
])

isoseqs_src_include_directories = include_directories('src/')

c3s_third_party_lib_impl = static_library(
  'c3s_third_party',
  c3s_third_party_cpp_sources,
  install : false,
  dependencies : c3s_lib_deps,
  include_directories : isoseqs_src_include_directories,
  cpp_args : c3s_flags)

pbsv_main = executable(
  'c3s', files([
    'src/main.cpp']),
  install : true,
  dependencies : c3s_lib_deps,
  link_with: c3s_third_party_lib_impl,
  cpp_args : c3s_flags)
