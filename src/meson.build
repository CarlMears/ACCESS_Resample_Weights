# ACCESS Resample Weights project
#
# Note that some suggested examples include:
#
# meson setup build_debug/ src/
# meson setup --buildtype release build_release/ src/
# meson setup -D b_sanitize=address,undefined build_debug_sanitize/ src/
# meson setup --buildtype release -D b_lto=true build_release_lto/ src/

project('Resample Weights for ACCESS',
        'fortran', 'c', 'cpp',
        version : '0.1',
        default_options : [
          'buildtype=debug',
          'warning_level=3',
          'b_pie=true',
          'c_std=c11',
          'fortran_std=f2008',
          'cpp_std=c++14',
        ],
        meson_version : '>=0.49')

if meson.get_compiler('fortran').get_argument_syntax() == 'gcc'
  add_project_arguments(
    # Always run the C pre-processor (for conditional compilation)
    '-cpp',
    '-ffpe-trap=zero,overflow',
    '-fno-plt',
    '-Wno-extra',
    # The code itself should already satisfy these next two, but this is
    # just to check that I didn't miss any
    '-fmodule-private',
    '-fimplicit-none',
    language : 'fortran')

  if meson.version().version_compare('<0.53')
    # This is instead handled by the "fortran_std" option for newer
    # versions of meson
    add_project_arguments('-std=f2008', language : 'fortran')
  endif

  add_project_arguments(
    '-Wformat=2',
    # There are some attributes used that only are supported in GCC
    # 10, so don't warn for unknown attributes
    '-Wno-attributes',
    language : 'c')

  add_project_link_arguments(
    '-Wl,-z,now',
    '-Wl,-z,relro',
    language : ['fortran', 'c'])

  if get_option('debug')
    add_project_arguments('-fcheck=all', language : 'fortran')
  endif
  if get_option('optimization') == '3'
    add_project_arguments('-ftree-vectorize', language : 'fortran')
  endif

elif meson.get_compiler('fortran').get_argument_syntax() == 'msvc'
  add_project_arguments(
    '/fpscomp:none',
    '/standard-semantics',
    '/assume:byterecl',
    '/fltconsistency',
    '/warn:all',
    '/check:all',
    '/traceback',
    '/Qip',
    '/Qipo',
    '/QxAVX',
    '/libs:dll',
    '/threads',
    language : 'fortran')
else
  warning('Not adding any compiler options')
endif

# -----------------------
# The netCDF library is used for reading and writing files. Both the C
# library and the Fortran bindings are required. The HDF5 library is
# indirectly used (though netCDF) but we need to do a version test to
# see if we need to add in calls to flock(2) or not.
if meson.version().version_compare('>=0.50')
  zlib = dependency('zlib', language : 'c')
  netcdf = dependency('netcdf', language : 'c')
  netcdff = dependency('netcdf', language : 'fortran')
  hdf5 = dependency('hdf5', language : 'c')
  hdf5_version = hdf5.version()
else
  zlib = dependency('zlib')
  netcdf = dependency('netcdf')
  netcdff = dependency('netcdf-fortran')
  hdf5 = dependency('hdf5', required : false)
  if hdf5.found()
    hdf5_version = hdf5.version()
  else
    # The EPEL 8 package doesn't include a pkg-config file and the
    # version of meson on CentOS 8 isn't new enough to be able to figure
    # out how to detect HDF5 without it. So try to grep the version from a file.
    message('Trying to detect HDF5 version')
    grep = find_program('grep')
    hdf5_version = run_command(grep, ['HDF5 Version', '/usr/lib64/libhdf5.settings']).stdout().strip().split(':')[1].strip()
  endif
endif

# -----------------------
# Generate the version information
git = find_program('git', required : false)
vcs_info = configuration_data()
if git.found()
  git_rev = run_command(git, [
                         '-C', meson.source_root(),
                         'log',
                         '--pretty=format:%H',
                         '-n', '1'
                       ]).stdout().strip()
  git_date = run_command(git, [
                          '-C', meson.source_root(),
                          'log',
                          '--pretty=format:%ai',
                          '-n', '1'
                        ]).stdout().strip()
  git_branch = run_command(git, [
                            '-C', meson.source_root(),
                            'symbolic-ref',
                            '--short', 'HEAD'
                          ]).stdout().strip()
  vcs_info.set('GIT_REV', git_rev)
  vcs_info.set('GIT_BRANCH', git_branch)
  vcs_info.set('GIT_DATE', git_date)
else
  vcs_info.set('GIT_REV', meson.project_version())
  vcs_info.set('GIT_BRANCH', 'unknown')
  vcs_info.set('GIT_DATE', 'unknown')
endif

log_lib = static_library('clogger', sources : ['syspgm/log.c'])
gzip = static_library('gzip', sources : ['syspgm/gzip_inflate.c','syspgm/zlib.f90'], dependencies: zlib)

executable('find_resampling_weights_precompute_g',
           sources : [
                'earth_grid/find_resampling_weights_precompute_g.f90',
                'earth_grid/amsr2_sensor_module.f90',
                'earth_grid/read_gobs_amsr2.f90',              
                'common/src/math/trig_degrees.f90',
                'common/src/simulation_param.f90',
                'common/src/resample/resampling_types.f90',
                'common/src/gain_patterns.f90',
                'common/src/geoloc/wgs84.f90',
                'common/src/math/svd_invert.f90',
                'common/src/math/pythag.f90'
           ],
           link_with : [gzip,log_lib],
           dependencies : [netcdff,zlib],
           install : true)

executable('find_resampling_weights_precompute_g_ssmi',
           sources : [
                'earth_grid/find_resampling_weights_precompute_g_ssmi.f90',
                'earth_grid/amsr2_sensor_module.f90',
                'earth_grid/ssmi_gains.f90',              
                'common/src/math/trig_degrees.f90',
                'common/src/simulation_param.f90',
                'common/src/resample/resampling_types.f90',
                'common/src/gain_patterns.f90',
                'common/src/geoloc/wgs84.f90',
                'common/src/math/svd_invert.f90',
                'common/src/math/pythag.f90'
           ],
           link_with : [gzip,log_lib],
           dependencies : [netcdff,zlib],
           install : true)

executable('find_resampling_weights_precompute_g_ssmis',
           sources : [
                'earth_grid/find_resampling_weights_precompute_g_ssmis.f90',
                'earth_grid/amsr2_sensor_module.f90',
                'earth_grid/ssmis_gains.f90',              
                'common/src/math/trig_degrees.f90',
                'common/src/simulation_param.f90',
                'common/src/resample/resampling_types.f90',
                'common/src/gain_patterns.f90',
                'common/src/geoloc/wgs84.f90',
                'common/src/math/svd_invert.f90',
                'common/src/math/pythag.f90'
           ],
           link_with : [gzip,log_lib],
           dependencies : [netcdff,zlib],
           install : true)

executable('calc_g_everywhere',
           sources : [
                'earth_grid/calc_g_everywhere.f90',
                'earth_grid/amsr2_sensor_module.f90',
                'earth_grid/ssmis_gains.f90', 
                'earth_grid/ssmi_gains.f90',
                'earth_grid/read_gobs_amsr2.f90',              
                'common/src/math/trig_degrees.f90',
                'common/src/simulation_param.f90',
                'common/src/resample/resampling_types.f90',
                'common/src/gain_patterns.f90',
                'common/src/geoloc/wgs84.f90',
                'common/src/math/svd_invert.f90',
                'common/src/math/pythag.f90'
           ],
           link_with : [gzip,log_lib],
           dependencies : [netcdff,zlib],
           install : true)
