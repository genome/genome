-- Revert remove_model_groups_from_analysis_project

BEGIN;

  ALTER TABLE config.analysis_project ADD COLUMN model_group_id character varying(64) NOT NULL;
  ALTER TABLE config.analysis_project ADD CONSTRAINT analysis_project_model_group FOREIGN KEY (model_group_id) REFERENCES model.model_group(id);

COMMIT;
