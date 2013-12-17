-- Deploy config.analysis_project.index_configuration_set_id
-- requires: analysis_project

BEGIN;

CREATE INDEX c_ap_configuration_set_id_index ON config.analysis_project USING btree (configuration_set_id);

COMMIT;
