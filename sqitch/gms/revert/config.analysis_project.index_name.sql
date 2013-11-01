-- Revert config.analysis_project.index_name

BEGIN;

DROP INDEX config.analysis_project_name_idx;

COMMIT;
