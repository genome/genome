-- Deploy config_analysis_menu_item_index
-- requires: config_analysis_menu_item

BEGIN;

CREATE INDEX analysis_menu_config_set_index ON config.analysis_menu_item USING btree (configuration_set_id);

COMMIT;
