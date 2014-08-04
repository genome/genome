-- Deploy config_tag_profile_item
-- requires: config_profile_item

BEGIN;
CREATE TABLE IF NOT EXISTS config.tag_profile_item (
    id character varying(64) PRIMARY KEY,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_by character varying(255) NOT NULL,
    tag_id character varying(64) NOT NULL REFERENCES config.tag(id),
    profile_item_id character varying(64) NOT NULL REFERENCES config.profile_item(id)
);

REVOKE ALL ON TABLE config.tag_profile_item FROM PUBLIC;
GRANT ALL ON TABLE config.tag_profile_item TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.tag_profile_item TO "gms-user";
COMMIT;
