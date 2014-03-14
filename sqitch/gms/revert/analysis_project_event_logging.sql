-- Revert analysis_project_event_logging

BEGIN;

  DROP TABLE IF EXISTS timeline.analysis_project;
  DROP TABLE IF EXISTS timeline.analysis_project_event_type;

COMMIT;
