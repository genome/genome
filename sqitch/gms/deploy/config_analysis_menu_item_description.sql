-- Deploy config_analysis_menu_item_description
-- requires: config_analysis_menu_item

BEGIN;

ALTER TABLE config.analysis_menu_item ADD COLUMN description VARCHAR(255);

COMMIT;
