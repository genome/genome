-- Deploy remove_unused_apconfig_tables
-- requires: config_set

BEGIN;
  ALTER TABLE config.analysis_project DROP COLUMN analysis_menu_item_id;
  ALTER TABLE config.analysis_project DROP COLUMN configuration_set_id;
  DROP TABLE config.analysis_menu_item;
  DROP TABLE config.set;
COMMIT;

