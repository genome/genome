-- Revert config.analysis_project.index_configuration_set_id

BEGIN;

DROP INDEX config.c_ap_configuration_set_id_index;

COMMIT;
