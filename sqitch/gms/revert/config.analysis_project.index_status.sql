-- Revert config.analysis_project.index_status

BEGIN;

DROP INDEX config.analysis_project_status_idx;

COMMIT;
