import logging
from galaxy.jobs.mapper import JobMappingException

log = logging.getLogger(__name__)

DESTINATION_IDS = {
    1 : 'vgl-all',
    2 : 'bigmem-all'
}
FAILURE_MESSAGE = 'This tool could not be run because of a misconfiguration in the Galaxy job running system, please report this error'


def dynamic_partition_cores(app, tool, job, user_email):
    destination = None
    destination_id = 'slurm'

    # build the param dictionary
    param_dict = job.get_param_values(app)

    if param_dict.get('__job_resource', {}).get('__job_resource__select') != 'yes':
        log.info("Job resource parameters not seleted, returning default destination")
        return destination_id

    # handle job resource parameters
    try:
        # validate params
        cores = int(param_dict['__job_resource']['cores'])
        partition = int(param_dict['__job_resource']['partition'])
        destination_id = DESTINATION_IDS[cores]
        destination = app.job_config.get_destination(destination_id)
        # set walltime
        if 'nativeSpecification' not in destination.params:
            destination.params['nativeSpecification'] = ''
        destination.params['nativeSpecification'] += ' --time=%s:00:00' % time
    except:
        # resource param selector not sent with tool form, job_conf.xml misconfigured
        log.warning('(%s) error, keys were: %s', job.id, param_dict.keys())
        raise JobMappingException(FAILURE_MESSAGE)

    log.info('returning destination: %s', destination_id)
    log.info('native specification: %s', destination.params.get('nativeSpecification'))
    return destination or destination_id
